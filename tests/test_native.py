"""Optional native integration tests. Set BLAST_BIN to a BLAST+ 2.17+ bin folder."""

import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


@unittest.skipUnless(os.environ.get("BLAST_BIN"), "Set BLAST_BIN to enable native megablast tests")
class NativeTests(unittest.TestCase):
    def test_corrupt_native_bundle_cannot_replace_existing_database(self):
        from apscale_blast2.db_home import install_precompiled_db_zip
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "db_previous"
            target.mkdir()
            sentinel = target / "KEEP"
            sentinel.write_text("previous database", encoding="ascii")
            archive = root / "invalid.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                for suffix in ["nin", "nhr", "nsq"]:
                    handle.writestr(f"db/db.{suffix}", "not a real index")
                handle.writestr("db_taxonomy.csv", "id,species\nref1,Testus alpha\n")
            executable = str(Path(os.environ["BLAST_BIN"]).resolve() / ("blastdbcmd.exe" if os.name == "nt" else "blastdbcmd"))
            with self.assertRaisesRegex(ValueError, "Cannot open"):
                install_precompiled_db_zip(str(archive), str(root), "previous", overwrite=True, blastdbcmd_exe=executable)
            self.assertEqual(sentinel.read_text(), "previous database")

    def test_build_and_cli_in_unicode_space_paths_without_stdin(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "project with spaces \u00f1"
            root.mkdir()
            randomizer = random.Random(812)
            sequence = "".join(randomizer.choice("ACGT") for _ in range(250))
            source = root / "reference.fa.gz"
            with gzip.open(source, "wt") as handle:
                handle.write(f">ref1;tax=k:Eukaryota,p:Phyluma,c:Classa,o:Ordera,f:Familya,g:Testus,s:Testus_alpha\n{sequence}\n")
            executable = lambda name: str(Path(os.environ["BLAST_BIN"]).resolve() / (name + (".exe" if os.name == "nt" else "")))
            env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONIOENCODING="cp1252", APSCALE_BLAST2_DB_HOME=str(root / "must_not_create"))
            command = [sys.executable, "-m", "apscale_blast2.cli"]
            build = command + ["build", "--recipe", "pr2", "--input", str(source), "--db-home", str(root / "dbs"), "--name", "test", "--makeblastdb-exe", executable("makeblastdb"), "--thresholds", "97,95,90,87,85"]
            result = subprocess.run(build, env=env, stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            bundle = root / "reference.zip"
            with zipfile.ZipFile(bundle, "w") as handle:
                for path in (root / "dbs" / "db_test").rglob("*"):
                    if path.is_file():
                        handle.write(path, path.relative_to(root / "dbs").as_posix())
            installed = subprocess.run(command + ["build", "--recipe", "precompiled", "--input", str(bundle), "--db-home", str(root / "copies"), "--name", "copy", "--blastdbcmd-exe", executable("blastdbcmd")], env=env, stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
            self.assertEqual(installed.returncode, 0, installed.stderr)
            query = root / "query.fasta.gz"
            with gzip.open(query, "wt") as handle:
                handle.write(f">001\n{sequence}\n>nohit\n{'N' * 250}\n")
            args = command + ["--fastas", str(query), "--db-for-all", str(root / "copies" / "db_copy"), "--out-dir", str(root / "results"), "--threads", "2", "--subset-size", "1", "--blastn-exe", executable("blastn"), "--makeblastdb-exe", executable("makeblastdb"), "--thresholds", "99,95,90,87,85", "--keep-tsv"]
            result = subprocess.run(args, env=env, stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            table = pd.read_excel(root / "results" / "taxonomy" / "query_taxonomy.xlsx", dtype={"unique ID": str}).fillna("")
            self.assertEqual(table["unique ID"].tolist(), ["001", "nohit"])
            self.assertEqual(table["Species"].tolist(), ["Testus alpha", "NoMatch"])
            info = json.loads((root / "results" / "taxonomy" / "query.runinfo.json").read_text())
            self.assertEqual(info["options"]["thresholds"], "99,95,90,87,85")
            self.assertEqual(info["options"]["task"], "megablast")
            self.assertTrue(list(Path(info["work_dir"]).rglob("*.tsv")))
            self.assertFalse((root / "must_not_create").exists())
            repeated = subprocess.run(args, env=env, stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn(b"--overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
