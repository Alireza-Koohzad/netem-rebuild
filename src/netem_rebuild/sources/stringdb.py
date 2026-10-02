from __future__ import annotations

import csv
import gzip
import hashlib
import tomllib
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATASET_CONFIG = PROJECT_ROOT / "configs" / "dataset_v1.toml"
SOURCE_CONFIG = PROJECT_ROOT / "configs" / "sources_v1.toml"

GENE_CATALOG = (
    PROJECT_ROOT / "data" / "interim" / "hgnc" / "gene_catalog.tsv"
)
CORE_GENES = (
    PROJECT_ROOT / "data" / "processed" / "nodes" / "genes.tsv"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "edges"
    / "gene_gene.tsv"
)


def load_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def verify_file(path: Path, expected_hash: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required source file not found: {path}")

    expected_hash = expected_hash.strip().lower()

    if len(expected_hash) != 64:
        raise ValueError(
            f"Invalid SHA-256 configured for {path.name}"
        )

    actual_hash = sha256(path)

    if actual_hash != expected_hash:
        raise ValueError(
            f"SHA-256 mismatch for {path.name}\n"
            f"Expected: {expected_hash}\n"
            f"Actual:   {actual_hash}"
        )


def load_gene_maps() -> tuple[
    set[str],
    set[str],
    dict[str, str],
    dict[str, str],
]:
    catalog_ids: set[str] = set()

    ensembl_candidates: dict[str, set[str]] = defaultdict(set)
    ncbi_candidates: dict[str, set[str]] = defaultdict(set)

    with GENE_CATALOG.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            hgnc_id = row["canonical_id"]
            catalog_ids.add(hgnc_id)

            ensembl_id = row["ensembl_gene_id"]
            if ensembl_id != "NA":
                ensembl_candidates[ensembl_id].add(hgnc_id)

            ncbi_id = row["ncbi_gene_id"]
            if ncbi_id != "NA":
                ncbi_candidates[ncbi_id].add(hgnc_id)

    # External IDs are only accepted when uniquely mapped
    # inside the frozen HGNC catalog.
    ensembl_to_hgnc = {
        external_id: next(iter(ids))
        for external_id, ids in ensembl_candidates.items()
        if len(ids) == 1
    }

    ncbi_to_hgnc = {
        external_id: next(iter(ids))
        for external_id, ids in ncbi_candidates.items()
        if len(ids) == 1
    }

    core_ids: set[str] = set()

    with CORE_GENES.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            core_ids.add(row["canonical_id"])

    return (
        catalog_ids,
        core_ids,
        ensembl_to_hgnc,
        ncbi_to_hgnc,
    )


def collect_physical_proteins(
    links_path: Path,
) -> set[str]:
    proteins: set[str] = set()

    with gzip.open(
        links_path,
        "rt",
        encoding="utf-8",
    ) as handle:
        header = next(handle).split()

        expected = {
            "protein1",
            "protein2",
            "experimental",
            "database",
            "textmining",
            "combined_score",
        }

        if set(header) != expected:
            raise ValueError(
                f"Unexpected STRING links columns: {header}"
            )

        for line in handle:
            protein1, protein2, *_ = line.split()
            proteins.update((protein1, protein2))

    return proteins


def build_protein_gene_mapping(
    aliases_path: Path,
    physical_proteins: set[str],
    catalog_ids: set[str],
    core_ids: set[str],
    ensembl_to_hgnc: dict[str, str],
    ncbi_to_hgnc: dict[str, str],
) -> tuple[
    dict[str, str],
    set[str],
    set[str],
    dict[str, str],
]:
    direct_hgnc: dict[str, set[str]] = defaultdict(set)
    ensembl_ids: dict[str, set[str]] = defaultdict(set)
    ncbi_ids: dict[str, set[str]] = defaultdict(set)

    with gzip.open(
        aliases_path,
        "rt",
        encoding="utf-8",
    ) as handle:
        header = next(handle).rstrip("\n").split("\t")

        expected = [
            "#string_protein_id",
            "alias",
            "source",
        ]

        if header != expected:
            raise ValueError(
                f"Unexpected STRING aliases columns: {header}"
            )

        for line in handle:
            protein, alias, source = (
                line.rstrip("\n").split("\t")
            )

            if protein not in physical_proteins:
                continue

            if source == "Ensembl_HGNC_hgnc_id":
                hgnc_id = (
                    alias
                    if alias.startswith("HGNC:")
                    else f"HGNC:{alias}"
                )
                direct_hgnc[protein].add(hgnc_id)

            elif source == "Ensembl_HGNC_ensembl_gene_id":
                ensembl_ids[protein].add(alias)

            elif source == "Ensembl_HGNC_entrez_id":
                ncbi_ids[protein].add(alias)

    resolved_core: dict[str, str] = {}
    ambiguous: set[str] = set()
    unmapped: set[str] = set()
    outside_core: dict[str, str] = {}

    for protein in physical_proteins:
        candidates: set[str] = set()

        candidates.update(
            hgnc_id
            for hgnc_id in direct_hgnc.get(protein, set())
            if hgnc_id in catalog_ids
        )

        candidates.update(
            ensembl_to_hgnc[ensembl_id]
            for ensembl_id in ensembl_ids.get(
                protein,
                set(),
            )
            if ensembl_id in ensembl_to_hgnc
        )

        candidates.update(
            ncbi_to_hgnc[ncbi_id]
            for ncbi_id in ncbi_ids.get(
                protein,
                set(),
            )
            if ncbi_id in ncbi_to_hgnc
        )

        if len(candidates) == 0:
            unmapped.add(protein)
            continue

        if len(candidates) > 1:
            ambiguous.add(protein)
            continue

        hgnc_id = next(iter(candidates))

        if hgnc_id in core_ids:
            resolved_core[protein] = hgnc_id
        else:
            outside_core[protein] = hgnc_id

    return (
        resolved_core,
        ambiguous,
        unmapped,
        outside_core,
    )


def build_gene_pairs(
    links_path: Path,
    protein_to_gene: dict[str, str],
) -> tuple[
    dict[tuple[str, str], dict[str, int]],
    dict[tuple[str, str], int],
    int,
    int,
    int,
]:
    best_evidence: dict[
        tuple[str, str],
        dict[str, int],
    ] = {}

    selected_protein_pair: dict[
        tuple[str, str],
        tuple[str, str],
    ] = {}

    support_count: dict[
        tuple[str, str],
        int,
    ] = defaultdict(int)

    raw_unique_protein_pairs = 0
    dropped_unresolved_pairs = 0
    collapsed_self_loops = 0

    with gzip.open(
        links_path,
        "rt",
        encoding="utf-8",
    ) as handle:
        header = next(handle).split()

        indices = {
            name: header.index(name)
            for name in (
                "experimental",
                "database",
                "textmining",
                "combined_score",
            )
        }

        for line in handle:
            parts = line.split()

            protein1 = parts[0]
            protein2 = parts[1]

            # The frozen STRING physical file was verified
            # to contain reciprocal A-B and B-A rows.
            # Keep one canonical orientation.
            if protein1 > protein2:
                continue

            raw_unique_protein_pairs += 1

            gene1 = protein_to_gene.get(protein1)
            gene2 = protein_to_gene.get(protein2)

            if gene1 is None or gene2 is None:
                dropped_unresolved_pairs += 1
                continue

            if gene1 == gene2:
                collapsed_self_loops += 1
                continue

            gene_pair = tuple(sorted((gene1, gene2)))
            protein_pair = (protein1, protein2)

            evidence = {
                "experimental_score": int(
                    parts[indices["experimental"]]
                ),
                "database_score": int(
                    parts[indices["database"]]
                ),
                "textmining_score": int(
                    parts[indices["textmining"]]
                ),
                "combined_score": int(
                    parts[indices["combined_score"]]
                ),
            }

            support_count[gene_pair] += 1

            previous = best_evidence.get(gene_pair)

            if previous is None:
                best_evidence[gene_pair] = evidence
                selected_protein_pair[gene_pair] = protein_pair
                continue

            if (
                evidence["combined_score"]
                > previous["combined_score"]
            ):
                best_evidence[gene_pair] = evidence
                selected_protein_pair[gene_pair] = protein_pair
                continue

            # Deterministic tie-break:
            # preserve evidence from one real protein pair.
            if (
                evidence["combined_score"]
                == previous["combined_score"]
                and protein_pair
                < selected_protein_pair[gene_pair]
            ):
                best_evidence[gene_pair] = evidence
                selected_protein_pair[gene_pair] = protein_pair

    return (
        best_evidence,
        support_count,
        raw_unique_protein_pairs,
        dropped_unresolved_pairs,
        collapsed_self_loops,
    )


def hgnc_sort_key(hgnc_id: str) -> int:
    return int(hgnc_id.split(":", 1)[1])


def main() -> None:
    dataset_config = load_toml(DATASET_CONFIG)
    source_config = load_toml(SOURCE_CONFIG)

    policy = dataset_config["gene_gene_network"]
    source = source_config["sources"]["string"]

    if policy["source"] != "STRING":
        raise ValueError(
            "Gene-gene network source must be STRING"
        )

    if policy["network_type"] != "physical":
        raise ValueError(
            "STRING network_type must be physical"
        )

    if policy["ambiguous_mapping_policy"] != "reject":
        raise ValueError(
            "Unsupported ambiguous mapping policy"
        )

    threshold = int(
        policy["primary_score_threshold"]
    )

    relation = policy["relation"]

    links_path = (
        PROJECT_ROOT / source["links_local_path"]
    )
    aliases_path = (
        PROJECT_ROOT / source["aliases_local_path"]
    )

    verify_file(
        links_path,
        source["links_sha256"],
    )
    verify_file(
        aliases_path,
        source["aliases_sha256"],
    )

    (
        catalog_ids,
        core_ids,
        ensembl_to_hgnc,
        ncbi_to_hgnc,
    ) = load_gene_maps()

    physical_proteins = collect_physical_proteins(
        links_path
    )

    (
        protein_to_gene,
        ambiguous,
        unmapped,
        outside_core,
    ) = build_protein_gene_mapping(
        aliases_path=aliases_path,
        physical_proteins=physical_proteins,
        catalog_ids=catalog_ids,
        core_ids=core_ids,
        ensembl_to_hgnc=ensembl_to_hgnc,
        ncbi_to_hgnc=ncbi_to_hgnc,
    )

    (
        best_evidence,
        support_count,
        raw_unique_protein_pairs,
        dropped_unresolved_pairs,
        collapsed_self_loops,
    ) = build_gene_pairs(
        links_path=links_path,
        protein_to_gene=protein_to_gene,
    )

    output_rows: list[dict[str, str | int]] = []

    for gene_pair, evidence in best_evidence.items():
        if evidence["combined_score"] < threshold:
            continue

        source_id, target_id = gene_pair

        output_rows.append(
            {
                "source_id": source_id,
                "target_id": target_id,
                "relation": relation,
                "combined_score": evidence[
                    "combined_score"
                ],
                "experimental_score": evidence[
                    "experimental_score"
                ],
                "database_score": evidence[
                    "database_score"
                ],
                "textmining_score": evidence[
                    "textmining_score"
                ],
                "supporting_protein_pairs": support_count[
                    gene_pair
                ],
            }
        )

    output_rows.sort(
        key=lambda row: (
            hgnc_sort_key(str(row["source_id"])),
            hgnc_sort_key(str(row["target_id"])),
        )
    )

    # --------------------------------------------------
    # Regression checks against the frozen V1 profile
    # --------------------------------------------------

    expected = policy["expected_build"]

    if (
        raw_unique_protein_pairs
        != expected["raw_unique_protein_pairs"]
    ):
        raise ValueError(
            "Unexpected raw unique STRING protein-pair count: "
            f"{raw_unique_protein_pairs}"
        )

    if (
        len(best_evidence)
        != expected["gene_pairs_before_threshold"]
    ):
        raise ValueError(
            "Unexpected pre-threshold gene-pair count: "
            f"{len(best_evidence)}"
        )

    if len(ambiguous) != expected["ambiguous_proteins"]:
        raise ValueError(
            "Unexpected ambiguous STRING protein count: "
            f"{len(ambiguous)}"
        )

    if len(output_rows) != expected["primary_edges"]:
        raise ValueError(
            "Unexpected final gene-gene edge count: "
            f"{len(output_rows)}"
        )

    active_genes = {
        str(row["source_id"])
        for row in output_rows
    } | {
        str(row["target_id"])
        for row in output_rows
    }

    if (
        len(active_genes)
        != expected["primary_active_genes"]
    ):
        raise ValueError(
            "Unexpected active gene count: "
            f"{len(active_genes)}"
        )

    # --------------------------------------------------
    # Final invariants
    # --------------------------------------------------

    edge_keys = {
        (
            row["source_id"],
            row["target_id"],
            row["relation"],
        )
        for row in output_rows
    }

    if len(edge_keys) != len(output_rows):
        raise ValueError(
            "Duplicate final gene-gene edges detected"
        )

    if any(
        row["source_id"] == row["target_id"]
        for row in output_rows
    ):
        raise ValueError(
            "Self-loop detected in final gene-gene layer"
        )

    if any(
        row["source_id"] not in core_ids
        or row["target_id"] not in core_ids
        for row in output_rows
    ):
        raise ValueError(
            "Final edge references a non-Core gene"
        )

    # --------------------------------------------------
    # Write output
    # --------------------------------------------------

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "source_id",
        "target_id",
        "relation",
        "combined_score",
        "experimental_score",
        "database_score",
        "textmining_score",
        "supporting_protein_pairs",
    ]

    with OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"STRING links verified: {links_path.name}")
    print(f"STRING aliases verified: {aliases_path.name}")
    print()

    print(
        "Physical STRING proteins:",
        f"{len(physical_proteins):,}",
    )
    print(
        "Resolved into Core V1:",
        f"{len(protein_to_gene):,}",
    )
    print(
        "Ambiguous mappings rejected:",
        f"{len(ambiguous):,}",
    )
    print(
        "Unmapped proteins:",
        f"{len(unmapped):,}",
    )
    print(
        "Mapped outside Core V1:",
        f"{len(outside_core):,}",
    )
    print()

    print(
        "Raw unique protein pairs:",
        f"{raw_unique_protein_pairs:,}",
    )
    print(
        "Dropped unresolved pairs:",
        f"{dropped_unresolved_pairs:,}",
    )
    print(
        "Collapsed self-loops:",
        f"{collapsed_self_loops:,}",
    )
    print(
        "Gene pairs before threshold:",
        f"{len(best_evidence):,}",
    )
    print()

    print(
        f"Primary threshold: {threshold}"
    )
    print(
        "Final Gene-Gene edges:",
        f"{len(output_rows):,}",
    )
    print(
        "Genes with >=1 Gene-Gene edge:",
        f"{len(active_genes):,}",
    )
    print(
        "Core gene coverage:",
        f"{100 * len(active_genes) / len(core_ids):.2f}%",
    )
    print()

    print(
        "Gene-Gene layer written:",
        OUTPUT.relative_to(PROJECT_ROOT),
    )


if __name__ == "__main__":
    main()