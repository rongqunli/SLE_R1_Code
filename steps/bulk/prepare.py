import csv
import hashlib
import os
from pathlib import Path
from configuration import TABLES

root = Path(os.environ["SLE_RUN_ROOT"]).resolve()
for relative, rows in TABLES.items():
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
for relative in ["SLE_v1/pipeline/3.2/data/annotation_library", "SLE_v1/pipeline/3.4/R_library"]:
    (root / relative).mkdir(parents=True, exist_ok=True)
raw = root / "SLE_v1/pipeline/3.1/data/raw"
path = root / "SLE_v1/pipeline/3.1/data/metadata/download_manifest.tsv"
path.parent.mkdir(parents=True, exist_ok=True)
with path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["file", "size_bytes", "sha256"])
    for name in ["GPL10558.annot.gz", "GSE65391_series_matrix.txt.gz"]:
        source = raw / name
        digest = hashlib.sha256()
        with source.open("rb") as raw_handle:
            for block in iter(lambda: raw_handle.read(1024 * 1024), b""):
                digest.update(block)
        writer.writerow([name, source.stat().st_size, digest.hexdigest()])
