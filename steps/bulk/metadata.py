import os
import shutil
import sys
from pathlib import Path

root = Path(os.environ["SLE_RUN_ROOT"]).resolve()
sys.path.insert(0, str(root / "SLE/metadata/scripts"))
import build_frozen_meta as build
import extract_geo_sample_metadata as extract

datasets = ["GSE65391", "GSE88884", "GSE88885", "GSE88886", "GSE135779", "GSE121239"]
build.INTERIM_DIR.mkdir(parents=True, exist_ok=True)
build.BY_DATASET_DIR.mkdir(parents=True, exist_ok=True)
all_records = []
for dataset in datasets:
    header = build.SOURCE_DIR / "geo_matrix_headers" / f"{dataset}_series_matrix.txt.header.txt"
    rows, summary = extract.parse_header(header)
    extract.write_tsv(build.INTERIM_DIR / f"{dataset}_series_matrix.txt.sample_metadata.tsv", rows)
    if dataset == "GSE135779":
        records, exclusions, supplement = build.gse135779_records()
        fields = list(dict.fromkeys((key for row in supplement for key in row)))
        build.write_tsv(
            build.OUTPUT_DIR / "GSE135779_published_supplement_extracted.tsv", supplement, fields
        )
    elif dataset.startswith("GSE8888"):
        records = build.gse8888x_records(dataset)
    else:
        records = getattr(build, dataset.lower() + "_records")()
    if len({record["sample_id"] for record in records}) != len(records):
        raise ValueError(f"Duplicate metadata sample IDs: {dataset}")
    build.write_tsv(build.BY_DATASET_DIR / f"{dataset}_meta.tsv", records, build.MASTER_FIELDS)
    all_records.extend(records)
build.write_tsv(build.OUTPUT_DIR / "master_sample_metadata.tsv", all_records, build.MASTER_FIELDS)
source = build.BY_DATASET_DIR / "GSE65391_meta.tsv"
target = root / "SLE_v1/pipeline/3.1/data/frozen/2026-07-25_v1.0/by_dataset/GSE65391_meta.tsv"
target.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(source, target)
for number in ["3.2", "3.3"]:
    target = root / f"SLE_v1/pipeline/{number}/data/metadata/GSE88884_meta.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(build.BY_DATASET_DIR / "GSE88884_meta.tsv", target)
