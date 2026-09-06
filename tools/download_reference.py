"""Download a public reference into ignored dbs/, recording provenance."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request
import urllib.parse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url")
    parser.add_argument("name")
    parser.add_argument("--md5")
    parser.add_argument("--sha256")
    parser.add_argument("--max-mib", type=int, default=512)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1] / "dbs" / "official"
    root.mkdir(parents=True, exist_ok=True)
    if Path(args.name).name != args.name or args.name in {".", ".."}:
        parser.error("name must be a filename")
    target = root / args.name
    part = target.with_name(target.name + ".partial")
    sha, md5, total = hashlib.sha256(), hashlib.md5(), 0
    request = urllib.request.Request(args.url, headers={"User-Agent": "apscale-blast2-reference-tests/1.0"})
    with urllib.request.urlopen(request, timeout=45) as response, part.open("wb") as handle:
        resolved_host = urllib.parse.urlsplit(response.url).netloc
        while block := response.read(1024 * 1024):
            total += len(block)
            if total > args.max_mib * 1024**2:
                raise ValueError("Download exceeds configured size limit")
            handle.write(block)
            sha.update(block)
            md5.update(block)
            if total % (16 * 1024**2) == 0:
                print(f"{args.name}: {total // 1024**2} MiB", flush=True)
    for expected, actual in [(args.md5, md5.hexdigest()), (args.sha256, sha.hexdigest())]:
        if expected and expected.lower() != actual:
            raise ValueError(f"Checksum mismatch for {args.name}")
    part.replace(target)
    metadata = {"url": args.url, "resolved_host": resolved_host, "downloaded_at": datetime.now(timezone.utc).isoformat(), "bytes": total, "sha256": sha.hexdigest(), "md5": md5.hexdigest(), "upstream_checksum_verified": bool(args.md5 or args.sha256)}
    target.with_name(target.name + ".provenance.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata), flush=True)


if __name__ == "__main__":
    main()
