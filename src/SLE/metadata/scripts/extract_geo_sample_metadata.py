from __future__ import annotations
import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def safe_name(value: str) -> str:
    value = value.strip().lower()
    value = value.replace("%", " percent ")
    value = re.sub("[^a-z0-9]+", "_", value)
    return value.strip("_") or "unnamed"


def deduplicate_columns(names: list[str]) -> list[str]:
    counts: dict[str, int] = defaultdict(int)
    result: list[str] = []
    for name in names:
        counts[name] += 1
        suffix = "" if counts[name] == 1 else f"__{counts[name]}"
        result.append(f"{name}{suffix}")
    return result


def parse_header(path: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    sample_rows: list[tuple[str, list[str]]] = []
    series_accession = ""
    platform_accessions: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!Series_geo_accession"):
                parsed = next(csv.reader([line], delimiter="\t", quotechar='"'))
                series_accession = parsed[1] if len(parsed) > 1 else ""
            elif line.startswith("!Series_platform_id"):
                parsed = next(csv.reader([line], delimiter="\t", quotechar='"'))
                platform_accessions.extend(parsed[1:])
            elif line.startswith("!Sample_"):
                parsed = next(csv.reader([line], delimiter="\t", quotechar='"'))
                sample_rows.append((parsed[0][len("!Sample_") :], parsed[1:]))
    accession_row = next((values for key, values in sample_rows if key == "geo_accession"), None)
    if not accession_row:
        raise ValueError(f"No !Sample_geo_accession row found in {path}")
    sample_count = len(accession_row)
    records: list[dict[str, str]] = [
        {
            "dataset": series_accession,
            "sample_id": gsm,
            "source_header_file": path.name,
            "platform_series": " | ".join(dict.fromkeys(platform_accessions)),
        }
        for gsm in accession_row
    ]
    characteristic_occurrences: dict[str, int] = defaultdict(int)
    ordinary_occurrences: dict[str, int] = defaultdict(int)
    for key, values in sample_rows:
        if len(values) != sample_count:
            raise ValueError(
                f"{path.name}: {key} has {len(values)} values; expected {sample_count}"
            )
        if key == "geo_accession":
            continue
        if key == "characteristics_ch1":
            parsed_values: list[tuple[str, str, str]] = []
            names: list[str] = []
            for raw in values:
                if ":" in raw:
                    label, value = raw.split(":", 1)
                    column = f"characteristic_{safe_name(label)}"
                    parsed_values.append((column, value.strip(), raw))
                else:
                    column = "characteristic_unlabelled"
                    parsed_values.append((column, raw.strip(), raw))
                names.append(column)
            if len(set(names)) == 1:
                base = names[0]
                characteristic_occurrences[base] += 1
                suffix = (
                    ""
                    if characteristic_occurrences[base] == 1
                    else f"__{characteristic_occurrences[base]}"
                )
                column = f"{base}{suffix}"
                for record, (_, value, _) in zip(records, parsed_values):
                    record[column] = value
            else:
                for record, (base, value, _) in zip(records, parsed_values):
                    if base not in record:
                        record[base] = value
                    else:
                        record[base] = f"{record[base]} | {value}"
            continue
        column_base = f"sample_{safe_name(key)}"
        ordinary_occurrences[column_base] += 1
        suffix = (
            ""
            if ordinary_occurrences[column_base] == 1
            else f"__{ordinary_occurrences[column_base]}"
        )
        column = f"{column_base}{suffix}"
        for record, value in zip(records, values):
            record[column] = value
    field_order = list(records[0])
    for record in records[1:]:
        for field in record:
            if field not in field_order:
                field_order.append(field)
    summary = {
        "dataset": series_accession,
        "source_header_file": path.name,
        "sample_count": sample_count,
        "platform_series": list(dict.fromkeys(platform_accessions)),
        "columns": field_order,
    }
    return (records, summary)


def write_tsv(path: Path, records: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for record in records:
        for key in record:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    all_records: list[dict[str, str]] = []
    for path in sorted(args.input_dir.glob("*.header.txt")):
        records, summary = parse_header(path)
        output_name = path.name.replace(".header.txt", ".sample_metadata.tsv")
        write_tsv(args.output_dir / output_name, records)
        summaries.append(summary)
        all_records.extend(records)
    write_tsv(args.output_dir / "all_geo_sample_metadata.tsv", all_records)
    (args.output_dir / "extraction_summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
