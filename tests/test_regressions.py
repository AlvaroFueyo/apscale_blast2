"""Fast regression suite; synthetic fixtures require no downloads or BLAST."""

import gzip
import itertools
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
import zipfile

import pandas as pd
import pyarrow.parquet as pq
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from apscale_blast2 import blast_runner as br
from apscale_blast2.assignment import TAX_COLS, assign_query
from apscale_blast2.cli import build_parser, discover_fastas, main, validate_range
from apscale_blast2.db_build_common import resolve_input_file, install_built_database
from apscale_blast2.db_build_diatbarcode import read_diatbarcode
from apscale_blast2.db_build_midori2 import midori2_taxonomy_table, _midori_token_to_name
from apscale_blast2.db_build_pr2 import pr2_taxonomy_table
from apscale_blast2.db_build_silva import silva_taxonomy_from_headers
from apscale_blast2.db_build_trnl import _load_taxonomy_table
from apscale_blast2.db_build_unite import unite_taxonomy_table
from apscale_blast2.db_defaults import defaults_path, save_db_defaults, get_thresholds_for_db
from apscale_blast2.db_home import install_precompiled_db_zip
from apscale_blast2.dbs import DatabaseSpec, ensure_db_prefix
from apscale_blast2.filtering import thresholds_to_dict
from apscale_blast2.io_utils import read_fasta_order, split_fasta, read_db_map
from apscale_blast2.streaming import TableSink, TAX_COLUMNS, run, _publish
from apscale_blast2.taxmap import load_taxmap_as_dict, find_taxmaps_paths
from apscale_blast2.taxonomy_clean import clean_species


def hit(subject="ref1", **changes):
    row = dict(zip(TAX_COLS, ["Eukaryota", "PhylumA", "ClassA", "OrderA", "FamilyA", "Testus", "Testus alpha"]))
    row.update({"Subject ID": subject, "Similarity": 100.0, "query_coverage": 100.0,
                "evalue": 1e-50, "taxonomy_mapped": True, "species_uncertainty": ""})
    row.update(changes)
    return row


class AssignmentTests(unittest.TestCase):
    def assign(self, rows, **options):
        return assign_query("query", rows, br.RunOptions(**options))

    def test_single_known_reference(self):
        result = self.assign([hit()])
        self.assertEqual((result["Species"], result["assigned_rank"], result["Flag"]), ("Testus alpha", "Species", ""))

    def test_no_match_has_explicit_status(self):
        self.assertEqual(self.assign([])["assignment_status"], "no_match")

    def test_missing_best_species_preserves_compatible_consensus(self):
        for rows in [[hit(Species=""), hit("ref2")], [hit(), hit("ref2", Species="")]]:
            self.assertEqual(self.assign(rows)["Species"], "Testus alpha")

    def test_all_blank_rank_does_not_hide_lower_conflict(self):
        result = self.assign([hit(Family=""), hit("ref2", Family="", Genus="Other", Species="Other beta")])
        self.assertTrue(result["Flag"].startswith("Fl2"))
        self.assertEqual(result["assigned_rank"], "Order")

    def test_all_blank_rank_does_not_truncate_compatible_lineage(self):
        self.assertEqual(self.assign([hit(Family=""), hit("ref2", Family="")])["Species"], "Testus alpha")

    def test_conflict_at_each_rank_truncates_below_it(self):
        for index, rank in enumerate(TAX_COLS[:-1]):
            with self.subTest(rank=rank):
                result = self.assign([hit(), hit("ref2", **{rank: "Other"})])
                self.assertTrue(result["Flag"].startswith(f"Fl{7-index}"))
                self.assertTrue(all(result[r] == "" for r in TAX_COLS[index:]))

    def test_two_species_label_is_not_a_species_assignment(self):
        result = self.assign([hit(), hit("ref2", Species="Testus beta")])
        self.assertEqual(result["Species"], "Testus alpha/beta")
        self.assertEqual(result["assigned_rank"], "Genus")
        self.assertEqual(result["assignment_status"], "ambiguous")

    def test_three_species_label(self):
        result = self.assign([hit(), hit("ref2", Species="Testus beta"), hit("ref3", Species="Testus gamma")])
        self.assertEqual(result["Species"], "Testus sp.")

    def test_permutation_invariance(self):
        rows = [hit(), hit("ref2", Species=""), hit("ref3", Species="Testus beta")]
        reference = self.assign(rows)
        for order in itertools.permutations(rows):
            self.assertEqual(self.assign(list(order)), reference)

    def test_soft_coverage_fallback(self):
        self.assertEqual(self.assign([hit(query_coverage=60.0)])["Species"], "Testus alpha")

    def test_soft_coverage_prefers_qualifying_hits(self):
        result = self.assign([hit(query_coverage=60.0), hit("ref2", Species="Testus beta", Similarity=99.0)])
        self.assertEqual(result["Species"], "Testus beta")

    def test_decimal_threshold_boundary(self):
        for identity, expected in [(99.5, "Testus alpha"), (99.49, "")]:
            self.assertEqual(self.assign([hit(Similarity=identity)], thresholds="99.5,95,90,87,85")["Species"], expected)

    def test_duplicate_reference_taxa_do_not_get_dominance(self):
        result = self.assign([hit(), hit("ref2"), hit("ref3", Species="Testus beta")])
        self.assertEqual(result["candidate_taxa"], 2)
        self.assertTrue(result["Flag"].startswith("Fl1"))

    def test_legacy_dominance(self):
        result = self.assign([hit(), hit("ref2"), hit("ref3", Species="Testus beta")], flag_scheme="apscale")
        self.assertTrue(result["Flag"].startswith("F1"))
        self.assertEqual(result["Species"], "Testus alpha")

    def test_missing_taxonomy_is_visible(self):
        self.assertEqual(self.assign([hit(taxonomy_mapped=False)])["assignment_status"], "taxonomy_missing")

    def test_qualified_reference_is_visible(self):
        result = self.assign([hit(Species="", species_uncertainty="qualified_name")])
        self.assertEqual(result["assigned_rank"], "Genus")
        self.assertEqual(result["reference_uncertainty"], "qualified_name")

    def test_hit_limit_is_visible(self):
        self.assertTrue(self.assign([hit()], max_target_seqs=1)["hit_limit_reached"])

    def test_uncertain_names_not_promoted_to_species(self):
        for name in ["Testus cf. alpha", "Testus aff. alpha", "Testus sp.", "Testus alpha x Testus beta", "Testus alpha/beta"]:
            with self.subTest(name=name):
                self.assertEqual(clean_species(name, genus="Testus"), "")

    def test_invalid_thresholds_rejected(self):
        for text in ["99,95", "99,,95,90,87,85", "NaN,95,90,87,85", "101,95,90,87,85", "90,95,90,87,85", "99,95,90,87,-1"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                thresholds_to_dict(text)

    def test_invalid_runtime_options_rejected(self):
        for kwargs in [dict(workers=0), dict(threads=-1), dict(subset_size=0), dict(max_target_seqs=1.5), dict(min_qcov=float("nan")), dict(filter_mode=3)]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                br.RunOptions(**kwargs)


class FileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def database(self, name="database", species="Testus alpha"):
        root = self.root / name
        for ext in ["nin", "nhr", "nsq"]:
            self.write(f"{name}/db/db.{ext}", "synthetic")
        self.write(f"{name}/db_taxonomy.csv", f"Accession,superkingdom,genus,species\nref1,Eukaryota,Testus,{species}\n")
        return root

    def test_gzip_fasta_and_single_file_discovery(self):
        path = self.root / "reads.fa.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(">a\nACGT\n>b\nACGU\n")
        self.assertEqual(read_fasta_order(str(path)), ["a", "b"])
        self.assertEqual(discover_fastas(str(path)), [str(path)])

    def test_invalid_fastas_fail_before_blast(self):
        for content in ["", "ACGT\n", ">a\n", ">a\nACGT\n>a\nACGT\n", ">a\nAXZ\n", ">\nACGT\n"]:
            with self.subTest(content=content), self.assertRaises(ValueError):
                read_fasta_order(str(self.write("bad.fa", content)))

    def test_split_no_empty_last_chunk(self):
        path = self.write("a.fa", ">a\nACGT\n>b\nACGT\n")
        chunks = split_fasta(str(path), str(self.root / "chunks"), 1)
        self.assertEqual([read_fasta_order(p) for p in chunks], [["a"], ["b"]])

    def test_local_taxonomy_never_uses_parent_parquet(self):
        database = self.database()
        pd.DataFrame({"id": ["ref1"], "species": ["Wrong beta"]}).to_parquet(self.root / "wrong_taxonomy.parquet")
        self.assertEqual(load_taxmap_as_dict(str(database / "db" / "db"))["ref1"][6], "Testus alpha")

    def test_missing_local_taxonomy_not_silently_borrowed(self):
        self.write("other_taxonomy.csv", "id,species\nref1,Testus alpha\n")
        directory = self.root / "unrelated"
        directory.mkdir()
        self.assertEqual(find_taxmaps_paths(str(directory)), [])

    def test_multiple_taxonomy_tables_rejected(self):
        self.write("a_taxonomy.csv", "id,species\nref1,Testus alpha\n")
        self.write("b_taxonomy.csv", "id,species\nref1,Testus alpha\n")
        with self.assertRaises(ValueError):
            find_taxmaps_paths(str(self.root))

    def test_numeric_identifier_leading_zeros_preserved(self):
        self.write("db_taxonomy.csv", "id,species\n001,Testus alpha\n")
        self.assertIn("001", load_taxmap_as_dict(str(self.root)))

    def test_conflicting_duplicate_ids_rejected(self):
        self.write("db_taxonomy.csv", "id,species\na,Testus alpha\na,Testus beta\n")
        with self.assertRaises(ValueError):
            load_taxmap_as_dict(str(self.root))

    def test_contradictory_species_and_genus_are_not_published_as_consistent(self):
        self.write("db_taxonomy.csv", "id,genus,species\na,Testus,Other alpha\n")
        row = load_taxmap_as_dict(str(self.root))["a"]
        self.assertEqual(row[5:8], ["Testus", "", "Other alpha"])
        self.assertIn("species/genus mismatch", row[9])

    def test_uncertainty_and_original_name_preserved(self):
        self.write("db_taxonomy.csv", "id,genus,species\na,Testus,Testus cf. alpha\nb,Testus,Testus alpha x Other beta\n")
        mapping = load_taxmap_as_dict(str(self.root))
        self.assertEqual(mapping["a"][6:8], ["", "Testus cf. alpha"])
        self.assertEqual(mapping["b"][5:7], ["", ""])
        self.assertEqual(mapping["b"][8], "hybrid")

    def test_database_prefix_exact_not_arbitrary_sibling(self):
        database = self.database()
        with self.assertRaises(ValueError):
            ensure_db_prefix(str(database / "db" / "wrong"))

    def test_explicit_prefix_with_sibling_alias(self):
        database = self.database()
        self.write("database/db/another.nal", "TITLE synthetic\n")
        self.assertEqual(ensure_db_prefix(str(database / "db" / "db")), str(database / "db" / "db"))

    def test_defaults_folder_and_prefix_use_same_location(self):
        database = self.database()
        self.assertEqual(defaults_path(str(database)), defaults_path(str(database / "db" / "db")))
        self.assertEqual(defaults_path(str(self.root / "future")), self.root / "future" / "apscale_blast2_defaults.json")

    def test_explicit_thresholds_override_database(self):
        database = self.database()
        save_db_defaults(str(database), "97,95,90,87,85")
        value, source, _ = get_thresholds_for_db(str(database), "99,95,90,87,85")
        self.assertEqual((value, source), ("99,95,90,87,85", "argument"))

    def test_taxonomy_cache_invalidated(self):
        database = self.database()
        logger = logging.getLogger("test")
        br.get_tax_dict_cached(str(database), logger)
        self.write("database/db_taxonomy.csv", "id,species\nref1,Other longername\n")
        self.assertEqual(br.get_tax_dict_cached(str(database), logger)["ref1"][6], "Other longername")

    def test_db_map_relative_paths_and_duplicate_rejection(self):
        mapping = self.write("map.csv", "fasta,db\na.fa,relative/db\n")
        # Windows temp directories may use an 8.3 alias for the same location.
        self.assertEqual(Path(read_db_map(str(mapping))["a.fa"]).resolve(),
                         (self.root / "relative" / "db").resolve())
        mapping.write_text("fasta,db\na.fa,one\na.fa,two\n")
        with self.assertRaises(ValueError):
            read_db_map(str(mapping))

    def archive(self, members):
        archive = self.root / "input.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            for name, data in members.items():
                handle.writestr(name, data)
        return archive

    def test_zip_traversal_rejected_without_touching_sibling(self):
        sentinel = self.write("sentinel", "KEEP")
        archive = self.archive({"../escape.fasta": ">a\nACGT\n"})
        with self.assertRaises(ValueError):
            resolve_input_file(str(archive), str(self.root))
        self.assertEqual(sentinel.read_text(), "KEEP")

    def test_zip_nested_fasta_extracted_to_owned_name(self):
        archive = self.archive({"a:b/../bad.fasta": ">a\nACGT\n"})
        with self.assertRaises(ValueError):
            resolve_input_file(str(archive), str(self.root))
        archive = self.archive({"directory/a.fasta": ">a\nACGT\n"})
        result = resolve_input_file(str(archive), str(self.root))
        self.assertEqual(Path(result).resolve(), (self.root / "a.fasta").resolve())
        self.assertTrue(Path(result).samefile(self.root / "a.fasta"))
        self.assertEqual(Path(result).read_text(encoding="utf-8"), ">a\nACGT\n")

    def test_zip_multiple_fasta_requires_explicit_selection(self):
        archive = self.archive({"a.fa": ">a\nACGT\n", "b.fa": ">b\nACGT\n"})
        with self.assertRaises(ValueError):
            resolve_input_file(str(archive), str(self.root))

    def test_invalid_precompiled_overwrite_preserves_previous_database(self):
        archive = self.archive({"db/db.nin": "invalid"})
        sentinel = self.write("home/db_previous/KEEP", "previous database")
        with self.assertRaises(ValueError):
            install_precompiled_db_zip(str(archive), str(self.root / "home"), name="previous", overwrite=True)
        self.assertTrue(sentinel.exists())

    def test_install_into_existing_empty_directory_not_nested(self):
        source = self.database()
        target = self.root / "target"
        target.mkdir()
        install_built_database(str(source), str(target))
        self.assertTrue((target / "db" / "db.nin").exists())
        self.assertFalse((target / "db_build").exists())

    def test_crux_headerless_first_row_kept(self):
        path = self.write("trnl.txt", "001\tEukaryota;P;NA;O;F;Testus;Testus alpha\n002\tEukaryota;P;NA;O;F;Testus;Testus beta\n")
        table = _load_taxonomy_table(str(path))
        self.assertEqual(table["Accession"].tolist(), ["001", "002"])

    def test_trnl_named_comma_csv_supported(self):
        path = self.write("trnl.csv", "Accession,superkingdom,phylum,class,order,family,genus,species\n001,Eukaryota,P,C,O,F,Testus,Testus alpha\n")
        self.assertEqual(_load_taxonomy_table(str(path)).iloc[0, 0], "001")

    def test_midori_preserves_qualifiers(self):
        self.assertEqual(_midori_token_to_name("Testus_cf._alpha_123"), "Testus cf. alpha")

    def test_unknown_header_formats_fail_explicitly(self):
        source = self.write("unknown.fa", ">plain\nACGT\n")
        for reader in [midori2_taxonomy_table, pr2_taxonomy_table, unite_taxonomy_table]:
            with self.subTest(reader=reader.__name__), self.assertRaises(ValueError):
                reader(str(source))

    def test_silva_requires_explicit_rank_map(self):
        source = self.write("silva.fa", ">id Eukaryota;Clade;P;C;O;F;Testus;Testus alpha\nACGT\n")
        with self.assertRaises(ValueError):
            silva_taxonomy_from_headers(str(source))
        ranks = self.write("ranks.txt", "Eukaryota;\t1\tdomain\nEukaryota;Clade;P;\t2\tphylum\nEukaryota;Clade;P;C;O;F;Testus;\t3\tgenus\n")
        table = silva_taxonomy_from_headers(str(source), str(ranks))
        self.assertEqual(table.iloc[0]["phylum"], "P")
        self.assertEqual(table.iloc[0]["class"], "")

    def test_excel_splits_sheets_and_keeps_formula_like_ids_as_text(self):
        path = self.root / "table.xlsx"
        with patch.object(TableSink, "excel_rows", 3):
            sink = TableSink(path, ["unique ID"], "excel")
            sink.append([{"unique ID": value} for value in ["=1+1", "002", "third"]])
            sink.close()
        book = load_workbook(path, read_only=True)
        self.addCleanup(book.close)
        self.assertEqual(len(book.worksheets), 2)
        self.assertEqual(book.worksheets[0]["A2"].data_type, "s")
        self.assertEqual(book.worksheets[0]["A2"].value, "=1+1")

    def test_empty_parquet_has_stable_schema(self):
        path = self.root / "table.parquet"
        sink = TableSink(path, TAX_COLUMNS, "parquet")
        sink.append([])
        sink.close()
        self.assertEqual(pq.read_schema(path).field("Similarity").type.__str__(), "double")

    def test_publication_rollback(self):
        a, b = self.write("a", "old a"), self.write("b", "old b")
        work = self.root / "work"
        work.mkdir()
        new_a = self.write("work/new_a", "new a")
        with self.assertRaises(FileNotFoundError):
            _publish([(new_a, a), (work / "missing", b)], work)
        self.assertEqual((a.read_text(), b.read_text()), ("old a", "old b"))

    def fake_blast(self, exe, prefix, query, *args, **kwargs):
        out = Path(query + ".tsv")
        out.write_text("".join(f"{qid}\tref1\tref1\tref1\t100\t1e-50\t100\t100\t0\t0\n" for qid in read_fasta_order(query)), encoding="utf-8")
        return str(out)

    def test_streaming_stable_across_chunks_and_workers(self):
        database = self.database()
        source = self.write("input.fa", "".join(f">q{i}\nACGT\n" for i in range(11)))
        tables = []
        with patch.object(br, "_run_single", side_effect=self.fake_blast):
            for index, (size, workers) in enumerate([(1, 1), (3, 2), (100, 1)]):
                opts = br.RunOptions(subset_size=size, workers=workers, threads=2, output_dir=str(self.root / f"result{index}"), output_format="parquet")
                result = run(str(source), str(self.root), DatabaseSpec(str(database)), opts)
                tables.append(pd.read_parquet(result["taxonomy_output"]))
        for table in tables[1:]:
            pd.testing.assert_frame_equal(table, tables[0])
        self.assertTrue(source.exists())

    def test_failed_overwrite_keeps_existing_outputs_and_releases_lock(self):
        database = self.database()
        source = self.write("input.fa", ">q1\nACGT\n")
        sentinel = self.write("result/taxonomy/input_taxonomy.xlsx", "KEEP")
        with patch.object(br, "_run_single", side_effect=RuntimeError("BLAST failed")) as search:
            with self.assertRaisesRegex(RuntimeError, "BLAST failed"):
                run(str(source), str(self.root), DatabaseSpec(str(database)), br.RunOptions(output_dir=str(self.root / "result"), overwrite=True))
            self.assertEqual(search.call_count, 1)
        self.assertEqual(sentinel.read_text(), "KEEP")
        self.assertFalse(list((self.root / "result").glob("*.lock")))

    def test_overwrite_requires_explicit_opt_in_before_blast(self):
        source = self.write("input.fa", ">q1\nACGT\n")
        self.write("result/taxonomy/input_taxonomy.xlsx", "KEEP")
        with patch.object(br, "_run_single") as search, self.assertRaises(FileExistsError):
            run(str(source), None, DatabaseSpec(str(self.root)), br.RunOptions(output_dir=str(self.root / "result")))
        search.assert_not_called()

    def test_cli_noninteractive_never_calls_input(self):
        database = self.database()
        source = self.write("input.fa", ">q1\nACGT\n")
        args = ["--fastas", str(source), "--db-for-all", str(database), "--out-dir", str(self.root / "result"), "--thresholds", "99,95,90,87,85"]
        with patch("builtins.input", side_effect=AssertionError("unexpected prompt")), patch("apscale_blast2.cli.require_blast_217"), patch.object(br, "_run_single", side_effect=self.fake_blast):
            self.assertEqual(main(args), 0)

    def test_cli_help_and_version_work_with_cp1252(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONIOENCODING="cp1252")
        for option in ["--help", "--version"]:
            result = subprocess.run([sys.executable, "-m", "apscale_blast2.cli", option], env=env, stdin=subprocess.DEVNULL, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_parser_help_and_integer_validation(self):
        self.assertIn("--thresholds", build_parser().format_help())
        self.assertEqual(validate_range("threads", 2, 0, 8, integer=True), 2)
        with self.assertRaises(ValueError):
            validate_range("threads", 2.5, 0, 8, integer=True)

    def diat_workbook(self, nodes, species="Testus alpha", sequence="AC-GT"):
        path = self.root / "diat.xlsx"
        with pd.ExcelWriter(path) as writer:
            pd.DataFrame({"Sequence ID": ["001"], "Species": [species], "Sequence": [sequence]}).to_excel(writer, sheet_name="sequences_info", index=False)
            tree = pd.DataFrame(nodes, columns=["taxon name", "rank", "parent taxon name"])
            for sheet in ["taxo_RCM", "taxo_Kociolek"]:
                tree.to_excel(writer, sheet_name=sheet, index=False)
        return path

    def test_diat_conflicting_parent_consensus_and_audit(self):
        nodes = [("Testus alpha", "species", "Testus"), ("Testus", "genus", "FamilyA"), ("Testus", "genus", "FamilyB"),
                 ("FamilyA", "family", "OrderA"), ("FamilyB", "family", "OrderA"), ("OrderA", "order", "")]
        path = self.diat_workbook(nodes)
        sequences, taxonomy = read_diatbarcode(str(path), "Kociolek")
        row = taxonomy.iloc[0]
        self.assertEqual((row["order"], row["family"], row["genus"], row["species"]), ("OrderA", "", "", ""))
        self.assertIn("FamilyA / FamilyB", row["taxonomy_conflict"])
        self.assertEqual(row["species_original"], "Testus alpha")
        self.assertEqual(sequences.iloc[0]["Sequence"], "ACGT")
        self.assertEqual(sequences.iloc[0]["gaps_removed"], 1)

    def test_diat_missing_ancestor_retains_known_ranks_and_audits(self):
        path = self.diat_workbook([("Testus alpha", "species", "Testus")])
        _, table = read_diatbarcode(str(path))
        self.assertEqual(table.iloc[0]["taxonomy_missing_nodes"], "Testus")
        self.assertEqual(table.iloc[0]["species"], "Testus alpha")
        self.assertEqual(table.iloc[0]["genus"], "")

    def test_diat_cycle_rejected(self):
        path = self.diat_workbook([("Testus alpha", "species", "Testus"), ("Testus", "genus", "Testus alpha")])
        with self.assertRaisesRegex(ValueError, "Cycle"):
            read_diatbarcode(str(path))

    def test_diat_legacy_blank_species_keeps_higher_ranks(self):
        path = self.root / "legacy.xlsx"
        pd.DataFrame({"Sequence ID": ["001"], "Sequence": ["ACGT"], "Species": [""], "Genus": ["Testus"],
                      "Subkingdom (following Algaebase 2018)": ["Eukaryota"], "Phylum (following Algaebase 2018)": ["P"],
                      "Class (following Round, Crawford & Mann 1990)": ["C"], "Order (following Round, Crawford & Mann 1990)": ["O"],
                      "Family (following Round, Crawford & Mann 1990)": ["F"]}).to_excel(path, sheet_name="diatbarcode v12", index=False)
        _, table = read_diatbarcode(str(path))
        self.assertEqual(table.iloc[0]["genus"], "Testus")
        self.assertEqual(table.iloc[0]["Accession"], "001")
        with self.assertRaises(ValueError):
            read_diatbarcode(str(path), "Kociolek")

    def test_pr2_stable_identifier_without_taxonomy_suffix(self):
        path = self.write("pr2.fa", ">abc;tax=k:Eukaryota,g:Testus,s:Testus_sp.\nACGT\n")
        table = pr2_taxonomy_table(str(path))
        self.assertEqual(table.iloc[0]["Accession"], "abc")
        self.assertEqual(br._resolve_sequence_id({"sseqid": "abc;tax=k:Eukaryota,g:Testus,s:Testus_sp"}, {"abc": ["known"]})[0], "abc")

    def test_build_cli_validation_never_prompts(self):
        with patch("builtins.input", side_effect=AssertionError("unexpected prompt")), self.assertRaises(SystemExit):
            main(["build", "--recipe", "silva", "--input", "file.fa"])

    def test_all_unmapped_hits_prevent_publication(self):
        database = self.database()
        source = self.write("input.fa", ">q1\nACGT\n")
        opts = br.RunOptions(output_dir=str(self.root / "result"), output_format="parquet")
        with patch.object(br, "_run_single", side_effect=self.fake_blast), patch.object(br, "get_tax_dict_cached", return_value={"other": ["Eukaryota"] * 7}), self.assertRaisesRegex(ValueError, "None of the BLAST hits"):
            run(str(source), None, DatabaseSpec(str(database)), opts)
        self.assertFalse((self.root / "result" / "taxonomy").exists())
        self.assertFalse(list((self.root / "result").glob("*.lock")))

    def test_empty_first_blast_chunk_does_not_break_parquet(self):
        database = self.database()
        source = self.write("input.fa", ">empty\nACGT\n>present\nACGT\n")

        def search(*args, **kwargs):
            out = self.fake_blast(*args, **kwargs)
            if read_fasta_order(args[2]) == ["empty"]:
                Path(out).write_text("")
            return out

        opts = br.RunOptions(output_dir=str(self.root / "result"), output_format="parquet", subset_size=1)
        with patch.object(br, "_run_single", side_effect=search):
            result = run(str(source), None, DatabaseSpec(str(database)), opts)
        table = pd.read_parquet(result["taxonomy_output"])
        self.assertEqual(table["assignment_status"].tolist(), ["no_match", "assigned"])

    def test_cancelled_search_does_not_start_process(self):
        cancel = threading.Event()
        cancel.set()
        with patch.object(br.subprocess, "Popen") as process, self.assertRaises(InterruptedError):
            br._run_single("blastn", "db", "query", 1, "megablast", 30, True, False, 50, .001, 50, cancel)
        process.assert_not_called()

    def test_worker_failure_cancels_other_inflight_work(self):
        database = self.database()
        source = self.write("input.fa", ">first\nACGT\n>fail\nACGT\n>waiting\nACGT\n>never\nACGT\n")
        cancelled = threading.Event()

        def search(*args, **kwargs):
            identifier = read_fasta_order(args[2])[0]
            if identifier == "fail":
                raise RuntimeError("test worker failure")
            if identifier == "waiting":
                self.assertTrue(kwargs["cancel_event"].wait(5), "Runner failed to cancel pending BLAST")
                cancelled.set()
                raise InterruptedError("cancelled")
            return self.fake_blast(*args, **kwargs)

        opts = br.RunOptions(output_dir=str(self.root / "result"), subset_size=1, workers=2, threads=2)
        started = time.monotonic()
        with patch.object(br, "_run_single", side_effect=search), self.assertRaisesRegex(RuntimeError, "test worker failure"):
            run(str(source), None, DatabaseSpec(str(database)), opts)
        self.assertLess(time.monotonic() - started, 4)

    def test_table_close_failure_does_not_publish_or_leak_lock(self):
        database = self.database()
        source = self.write("input.fa", ">q1\nACGT\n")
        close = TableSink.close

        def failed_close(sink):
            close(sink)
            raise OSError("disk write failure")

        opts = br.RunOptions(output_dir=str(self.root / "result"), output_format="parquet")
        with patch.object(br, "_run_single", side_effect=self.fake_blast), patch.object(TableSink, "close", failed_close), self.assertRaisesRegex(OSError, "disk write failure"):
            run(str(source), None, DatabaseSpec(str(database)), opts)
        self.assertFalse((self.root / "result" / "taxonomy").exists())
        self.assertFalse(list((self.root / "result").glob("*.lock")))


if __name__ == "__main__":
    unittest.main()
