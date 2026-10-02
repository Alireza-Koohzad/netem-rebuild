from __future__ import annotations

import csv
import hashlib
import tomllib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATASET_CONFIG = PROJECT_ROOT / "configs" / "dataset_v1.toml"
SOURCE_CONFIG = PROJECT_ROOT / "configs" / "sources_v1.toml"

CATALOG_OUTPUT = PROJECT_ROOT / "data" / "interim" / "hgnc" / "gene_catalog.tsv"
CORE_OUTPUT = PROJECT_ROOT / "data" / "processed" / "nodes" / "genes.tsv"


def load_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def normalize(value: str | None) -> str:
    value = (value or "").strip()
    return value if value else "NA"


def main() -> None:
    dataset_config = load_toml(DATASET_CONFIG)
    source_config = load_toml(SOURCE_CONFIG)

    gene_policy = dataset_config["gene_universe"]
    hgnc_source = source_config["sources"]["hgnc"]

    raw_path = PROJECT_ROOT / hgnc_source["local_path"]

    if not raw_path.exists():
        raise FileNotFoundError(f"HGNC raw file not found: {raw_path}")

    expected_hash = hgnc_source["sha256"].strip().lower()

    if not expected_hash:
        raise ValueError("HGNC SHA-256 is missing in configs/sources_v1.toml")

    actual_hash = sha256(raw_path)

    if actual_hash != expected_hash:
        raise ValueError(
            "HGNC SHA-256 mismatch.\n"
            f"Expected: {expected_hash}\n"
            f"Actual:   {actual_hash}"
        )

    required_source_columns = {
        "hgnc_id",
        "symbol",
        "name",
        "locus_group",
        "locus_type",
        "status",
        "entrez_id",
        "ensembl_gene_id",
    }

    with raw_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        missing_columns = required_source_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                f"HGNC source is missing columns: {sorted(missing_columns)}"
            )

        rows = list(reader)

    catalog_rows: list[dict[str, str]] = []
    core_rows: list[dict[str, str]] = []

    seen_ids: set[str] = set()

    for row in rows:
        hgnc_id = row["hgnc_id"].strip()

        if not hgnc_id:
            raise ValueError("HGNC record with empty hgnc_id found")

        if hgnc_id in seen_ids:
            raise ValueError(f"Duplicate HGNC ID found: {hgnc_id}")

        seen_ids.add(hgnc_id)

        if row["status"].strip() != "Approved":
            continue

        record = {
            "canonical_id": hgnc_id,
            "symbol": row["symbol"].strip(),
            "name": row["name"].strip(),
            "locus_group": row["locus_group"].strip(),
            "locus_type": row["locus_type"].strip(),
            "ensembl_gene_id": normalize(row["ensembl_gene_id"]),
            "ncbi_gene_id": normalize(row["entrez_id"]),
        }

        catalog_rows.append(record)

        if (
            record["locus_group"] == gene_policy["required_locus_group"]
            and record["locus_type"] == gene_policy["required_locus_type"]
        ):
            core_rows.append(
                {
                    "canonical_id": record["canonical_id"],
                    "symbol": record["symbol"],
                    "name": record["name"],
                    "ensembl_gene_id": record["ensembl_gene_id"],
                    "ncbi_gene_id": record["ncbi_gene_id"],
                }
            )

    catalog_rows.sort(key=lambda row: int(row["canonical_id"].split(":")[1]))
    core_rows.sort(key=lambda row: int(row["canonical_id"].split(":")[1]))

    CATALOG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    CORE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    catalog_fields = [
        "canonical_id",
        "symbol",
        "name",
        "locus_group",
        "locus_type",
        "ensembl_gene_id",
        "ncbi_gene_id",
    ]

    with CATALOG_OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=catalog_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(catalog_rows)

    core_fields = [
        "canonical_id",
        "symbol",
        "name",
        "ensembl_gene_id",
        "ncbi_gene_id",
    ]

    with CORE_OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=core_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(core_rows)

    print(f"HGNC source verified: {raw_path.name}")
    print(f"Approved HGNC catalog: {len(catalog_rows):,}")
    print(f"Core V1 protein-coding genes: {len(core_rows):,}")
    print(f"Catalog written: {CATALOG_OUTPUT.relative_to(PROJECT_ROOT)}")
    print(f"Core genes written: {CORE_OUTPUT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()