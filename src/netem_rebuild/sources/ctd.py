from __future__ import annotations

import csv
import gzip
import hashlib
import tomllib
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATASET_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "dataset_v1.toml"
)

SOURCE_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "sources_v1.toml"
)

CORE_GENES = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "nodes"
    / "genes.tsv"
)

CHEMICALS_OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "nodes"
    / "chemicals.tsv"
)

EDGES_OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "edges"
    / "chemical_gene.tsv"
)


CTD_FIELDS = [
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "GeneForms",
    "Organism",
    "OrganismID",
    "Interaction",
    "InteractionActions",
    "PubMedIDs",
]


def load_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def verify_file(
    path: Path,
    expected_hash: str,
) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required CTD source file not found: {path}"
        )

    expected_hash = expected_hash.strip().lower()

    if len(expected_hash) != 64:
        raise ValueError(
            "Configured CTD SHA-256 must contain "
            "exactly 64 hexadecimal characters"
        )

    actual_hash = sha256(path)

    if actual_hash != expected_hash:
        raise ValueError(
            "CTD SHA-256 mismatch.\n"
            f"Expected: {expected_hash}\n"
            f"Actual:   {actual_hash}"
        )


def inspect_metadata(
    path: Path,
) -> tuple[str, list[str]]:
    report_created = ""
    fields: list[str] = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8-sig",
    ) as handle:
        for line in handle:
            stripped = line.strip()

            if stripped.startswith("# Report created:"):
                report_created = stripped.split(
                    ":",
                    1,
                )[1].strip()

            elif stripped.startswith("# ChemicalName,"):
                fields = stripped[2:].split(",")

            if report_created and fields:
                break

    if not report_created:
        raise ValueError(
            "CTD metadata does not contain "
            "'# Report created:'"
        )

    if fields != CTD_FIELDS:
        raise ValueError(
            "Unexpected CTD field definition.\n"
            f"Expected: {CTD_FIELDS}\n"
            f"Actual:   {fields}"
        )

    return report_created, fields


def load_core_gene_mapping() -> tuple[
    dict[str, str],
    set[str],
]:
    ncbi_candidates: dict[
        str,
        set[str],
    ] = defaultdict(set)

    core_ids: set[str] = set()

    if not CORE_GENES.exists():
        raise FileNotFoundError(
            f"Core gene file not found: {CORE_GENES}"
        )

    with CORE_GENES.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required_columns = {
            "canonical_id",
            "ncbi_gene_id",
        }

        missing_columns = (
            required_columns
            - set(reader.fieldnames or [])
        )

        if missing_columns:
            raise ValueError(
                "Core gene file is missing columns: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            hgnc_id = row[
                "canonical_id"
            ].strip()

            core_ids.add(hgnc_id)

            ncbi_id = row[
                "ncbi_gene_id"
            ].strip()

            if (
                ncbi_id
                and ncbi_id != "NA"
            ):
                ncbi_candidates[
                    ncbi_id
                ].add(hgnc_id)

    # Only uniquely mapped external identifiers
    # are accepted.
    ncbi_to_hgnc = {
        ncbi_id: next(iter(hgnc_ids))
        for ncbi_id, hgnc_ids
        in ncbi_candidates.items()
        if len(hgnc_ids) == 1
    }

    return (
        ncbi_to_hgnc,
        core_ids,
    )


def split_pipe_values(
    value: str,
) -> set[str]:
    return {
        item.strip()
        for item in value.split("|")
        if item.strip()
    }


def is_complex_interaction(
    interaction: str,
) -> bool:
    return (
        "[" in interaction
        or "]" in interaction
    )


def contains_cotreatment(
    interaction_actions: str,
) -> bool:
    return (
        "affects^cotreatment"
        in split_pipe_values(
            interaction_actions
        )
    )


def contains_no_effect(
    interaction: str,
    interaction_actions: str,
) -> bool:
    interaction_lower = (
        interaction.lower()
    )

    actions_lower = (
        interaction_actions.lower()
    )

    return (
        "does not affect"
        in interaction_lower
        or "does not affect"
        in actions_lower
    )


def mesh_id(
    raw_chemical_id: str,
) -> str:
    return f"MESH:{raw_chemical_id}"


def chemical_sort_key(
    canonical_id: str,
) -> tuple[str, str]:
    raw_id = canonical_id.split(
        ":",
        1,
    )[1]

    return (
        raw_id[:1],
        raw_id,
    )


def hgnc_sort_key(
    hgnc_id: str,
) -> int:
    return int(
        hgnc_id.split(
            ":",
            1,
        )[1]
    )


def pubmed_sort_key(
    value: str,
) -> tuple[int, int | str]:
    if value.isdigit():
        return (
            0,
            int(value),
        )

    return (
        1,
        value,
    )


def main() -> None:
    dataset_config = load_toml(
        DATASET_CONFIG
    )

    source_config = load_toml(
        SOURCE_CONFIG
    )

    policy = dataset_config[
        "chemical_gene_network"
    ]

    expected = policy[
        "expected_build"
    ]

    source = source_config[
        "sources"
    ]["ctd"]

    # --------------------------------------------------
    # Validate configured policy
    # --------------------------------------------------

    if policy["source"] != "CTD":
        raise ValueError(
            "Chemical-gene source must be CTD"
        )

    if (
        policy["chemical_namespace"]
        != "MESH"
    ):
        raise ValueError(
            "CTD chemical namespace must be MESH"
        )

    if (
        policy["gene_mapping"]
        != "NCBI_Gene_to_HGNC"
    ):
        raise ValueError(
            "Unsupported CTD gene mapping policy"
        )

    if not policy["binary_only"]:
        raise ValueError(
            "CTD V1 requires binary_only=true"
        )

    if not policy[
        "exclude_cotreatment"
    ]:
        raise ValueError(
            "CTD V1 requires "
            "exclude_cotreatment=true"
        )

    if not policy[
        "exclude_no_effect"
    ]:
        raise ValueError(
            "CTD V1 requires "
            "exclude_no_effect=true"
        )

    if not policy["directed"]:
        raise ValueError(
            "CTD V1 Chemical-Gene "
            "storage direction must be directed"
        )

    if (
        policy["canonical_direction"]
        != "chemical_to_gene"
    ):
        raise ValueError(
            "Unsupported CTD canonical direction"
        )

    organism_id = str(
        policy["organism_id"]
    )

    relation = policy[
        "relation"
    ]

    raw_path = (
        PROJECT_ROOT
        / source["local_path"]
    )

    # --------------------------------------------------
    # Verify raw source
    # --------------------------------------------------

    verify_file(
        raw_path,
        source["sha256"],
    )

    report_created, _ = (
        inspect_metadata(
            raw_path
        )
    )

    # Our frozen source config says
    # report_created = 2026-09-29.
    # Verify this against the human-readable
    # metadata embedded in the CTD file.
    expected_report_date = source[
        "report_created"
    ]

    if expected_report_date != "2026-09-29":
        raise ValueError(
            "Unexpected configured CTD "
            f"report date: {expected_report_date}"
        )

    if (
        "Sep 29" not in report_created
        or "2026" not in report_created
    ):
        raise ValueError(
            "Unexpected CTD report creation "
            f"metadata: {report_created}"
        )

    # --------------------------------------------------
    # Core gene mapping
    # --------------------------------------------------

    (
        ncbi_to_hgnc,
        core_ids,
    ) = load_core_gene_mapping()

    # --------------------------------------------------
    # Chemical identity
    # --------------------------------------------------

    chemical_names: dict[
        str,
        set[str],
    ] = defaultdict(set)

    # CAS RN is optional metadata.
    # A valid MeSH chemical concept does not
    # necessarily have a CAS registry number.
    chemical_cas: dict[
        str,
        set[str],
    ] = defaultdict(set)

    # --------------------------------------------------
    # Edge evidence
    # --------------------------------------------------

    curated_rows: dict[
        tuple[str, str],
        int,
    ] = defaultdict(int)

    pubmed_ids: dict[
        tuple[str, str],
        set[str],
    ] = defaultdict(set)

    interaction_actions: dict[
        tuple[str, str],
        set[str],
    ] = defaultdict(set)

    gene_forms: dict[
        tuple[str, str],
        set[str],
    ] = defaultdict(set)

    # --------------------------------------------------
    # Build counters
    # --------------------------------------------------

    human_rows = 0
    mapped_human_rows = 0
    strict_curated_rows = 0

    complex_rows_rejected = 0
    cotreatment_rows_rejected = 0
    no_effect_rows_rejected = 0
    unmapped_gene_rows = 0

    # --------------------------------------------------
    # Stream CTD
    # --------------------------------------------------

    with gzip.open(
        raw_path,
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        data_lines = (
            line
            for line in handle
            if (
                line.strip()
                and not line.startswith("#")
            )
        )

        reader = csv.reader(
            data_lines
        )

        for values in reader:
            if (
                len(values)
                != len(CTD_FIELDS)
            ):
                raise ValueError(
                    "Unexpected CTD data-row "
                    "column count: "
                    f"{len(values)}"
                )

            row = dict(
                zip(
                    CTD_FIELDS,
                    values,
                )
            )

            # ------------------------------------------
            # Human-only filter
            # ------------------------------------------

            if (
                row["OrganismID"].strip()
                != organism_id
            ):
                continue

            human_rows += 1

            # ------------------------------------------
            # Gene mapping
            # ------------------------------------------

            gene_id = row[
                "GeneID"
            ].strip()

            hgnc_id = (
                ncbi_to_hgnc.get(
                    gene_id
                )
            )

            if hgnc_id is None:
                unmapped_gene_rows += 1
                continue

            mapped_human_rows += 1

            # ------------------------------------------
            # Chemical identity
            # ------------------------------------------

            chemical_id = row[
                "ChemicalID"
            ].strip()

            if not chemical_id:
                raise ValueError(
                    "CTD human mapped row "
                    "contains empty ChemicalID"
                )

            chemical = mesh_id(
                chemical_id
            )

            # ------------------------------------------
            # Strict V1 evidence policy
            # ------------------------------------------

            interaction = row[
                "Interaction"
            ].strip()

            actions_text = row[
                "InteractionActions"
            ].strip()

            if (
                policy["binary_only"]
                and is_complex_interaction(
                    interaction
                )
            ):
                complex_rows_rejected += 1
                continue

            if (
                policy[
                    "exclude_cotreatment"
                ]
                and contains_cotreatment(
                    actions_text
                )
            ):
                cotreatment_rows_rejected += 1
                continue

            if (
                policy[
                    "exclude_no_effect"
                ]
                and contains_no_effect(
                    interaction,
                    actions_text,
                )
            ):
                no_effect_rows_rejected += 1
                continue

            strict_curated_rows += 1

            # ------------------------------------------
            # Chemical attributes
            # ------------------------------------------

            name = row[
                "ChemicalName"
            ].strip()

            cas_rn = row[
                "CasRN"
            ].strip()

            if not name:
                raise ValueError(
                    "Missing ChemicalName "
                    f"for {chemical}"
                )

            chemical_names[
                chemical
            ].add(name)

            # CAS is optional.
            if cas_rn:
                chemical_cas[
                    chemical
                ].add(cas_rn)

            # ------------------------------------------
            # Edge evidence
            # ------------------------------------------

            pair = (
                chemical,
                hgnc_id,
            )

            curated_rows[
                pair
            ] += 1

            pubmed_ids[
                pair
            ].update(
                split_pipe_values(
                    row["PubMedIDs"]
                )
            )

            interaction_actions[
                pair
            ].update(
                split_pipe_values(
                    actions_text
                )
            )

            gene_forms[
                pair
            ].update(
                split_pipe_values(
                    row["GeneForms"]
                )
            )

    # --------------------------------------------------
    # Validate chemical identity
    # --------------------------------------------------

    for (
        chemical,
        names,
    ) in chemical_names.items():
        if len(names) != 1:
            raise ValueError(
                "ChemicalID maps to multiple "
                "names: "
                f"{chemical} -> "
                f"{sorted(names)}"
            )

    for chemical in chemical_names:
        cas_values = (
            chemical_cas.get(
                chemical,
                set(),
            )
        )

        if len(cas_values) > 1:
            raise ValueError(
                "ChemicalID maps to multiple "
                "nonblank CAS values: "
                f"{chemical} -> "
                f"{sorted(cas_values)}"
            )

    # --------------------------------------------------
    # Build chemical nodes
    # --------------------------------------------------

    chemical_rows: list[
        dict[str, str]
    ] = []

    for chemical in chemical_names:
        names = chemical_names[
            chemical
        ]

        cas_values = (
            chemical_cas.get(
                chemical,
                set(),
            )
        )

        chemical_rows.append(
            {
                "canonical_id":
                    chemical,
                "name":
                    next(iter(names)),
                "cas_rn":
                    (
                        next(
                            iter(
                                cas_values
                            )
                        )
                        if cas_values
                        else "NA"
                    ),
            }
        )

    chemical_rows.sort(
        key=lambda row:
        chemical_sort_key(
            row["canonical_id"]
        )
    )

    missing_cas_chemicals = sum(
        row["cas_rn"] == "NA"
        for row in chemical_rows
    )

    # --------------------------------------------------
    # Build Chemical-Gene edges
    # --------------------------------------------------

    edge_rows: list[
        dict[str, str | int]
    ] = []

    for pair in curated_rows:
        chemical, hgnc_id = pair

        pmids = sorted(
            pubmed_ids[pair],
            key=pubmed_sort_key,
        )

        actions = sorted(
            interaction_actions[
                pair
            ]
        )

        forms = sorted(
            gene_forms[
                pair
            ]
        )

        edge_rows.append(
            {
                "source_id":
                    chemical,
                "target_id":
                    hgnc_id,
                "relation":
                    relation,
                "supporting_curated_rows":
                    curated_rows[pair],
                "supporting_pubmed_count":
                    len(pmids),
                "pubmed_ids":
                    "|".join(pmids),
                "interaction_actions":
                    (
                        "|".join(
                            actions
                        )
                        if actions
                        else "NA"
                    ),
                "gene_forms":
                    (
                        "|".join(
                            forms
                        )
                        if forms
                        else "NA"
                    ),
            }
        )

    edge_rows.sort(
        key=lambda row: (
            chemical_sort_key(
                str(
                    row[
                        "source_id"
                    ]
                )
            ),
            hgnc_sort_key(
                str(
                    row[
                        "target_id"
                    ]
                )
            ),
        )
    )

    active_core_genes = {
        str(
            row["target_id"]
        )
        for row in edge_rows
    }

    active_chemicals = {
        str(
            row["source_id"]
        )
        for row in edge_rows
    }

    # --------------------------------------------------
    # Regression checks
    # --------------------------------------------------

    if (
        strict_curated_rows
        != expected[
            "strict_curated_rows"
        ]
    ):
        raise ValueError(
            "Unexpected strict CTD "
            "row count: "
            f"{strict_curated_rows}"
        )

    if (
        len(chemical_rows)
        != expected["chemicals"]
    ):
        raise ValueError(
            "Unexpected CTD chemical "
            "count: "
            f"{len(chemical_rows)}"
        )

    if (
        len(edge_rows)
        != expected["edges"]
    ):
        raise ValueError(
            "Unexpected CTD edge "
            "count: "
            f"{len(edge_rows)}"
        )

    if (
        len(active_core_genes)
        != expected[
            "active_core_genes"
        ]
    ):
        raise ValueError(
            "Unexpected CTD active "
            "Core gene count: "
            f"{len(active_core_genes)}"
        )

    if (
        len(active_chemicals)
        != len(chemical_rows)
    ):
        raise ValueError(
            "At least one retained "
            "chemical has no final edge"
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
        for row in edge_rows
    }

    if len(edge_keys) != len(edge_rows):
        raise ValueError(
            "Duplicate final "
            "Chemical-Gene edges detected"
        )

    if any(
        row["target_id"]
        not in core_ids
        for row in edge_rows
    ):
        raise ValueError(
            "Chemical-Gene edge "
            "references a non-Core gene"
        )

    chemical_ids = {
        row["canonical_id"]
        for row in chemical_rows
    }

    if (
        len(chemical_ids)
        != len(chemical_rows)
    ):
        raise ValueError(
            "Duplicate final chemical "
            "node IDs detected"
        )

    if any(
        row["source_id"]
        not in chemical_ids
        for row in edge_rows
    ):
        raise ValueError(
            "Chemical-Gene edge "
            "references an unknown chemical"
        )

    if any(
        not str(
            row["source_id"]
        ).startswith("MESH:")
        for row in edge_rows
    ):
        raise ValueError(
            "Non-MESH chemical "
            "identifier detected"
        )

    if any(
        int(
            row[
                "supporting_curated_rows"
            ]
        ) < 1
        for row in edge_rows
    ):
        raise ValueError(
            "Invalid curated-row "
            "support count"
        )

    if any(
        int(
            row[
                "supporting_pubmed_count"
            ]
        ) < 1
        for row in edge_rows
    ):
        raise ValueError(
            "Chemical-Gene edge "
            "without PMID evidence"
        )

    # --------------------------------------------------
    # Write chemical nodes
    # --------------------------------------------------

    CHEMICALS_OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    chemical_fields = [
        "canonical_id",
        "name",
        "cas_rn",
    ]

    with CHEMICALS_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=chemical_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerows(
            chemical_rows
        )

    # --------------------------------------------------
    # Write Chemical-Gene edges
    # --------------------------------------------------

    EDGES_OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    edge_fields = [
        "source_id",
        "target_id",
        "relation",
        "supporting_curated_rows",
        "supporting_pubmed_count",
        "pubmed_ids",
        "interaction_actions",
        "gene_forms",
    ]

    with EDGES_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=edge_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerows(
            edge_rows
        )

    # --------------------------------------------------
    # Build summary
    # --------------------------------------------------

    coverage = (
        100
        * len(active_core_genes)
        / len(core_ids)
    )

    print(
        "CTD source verified:",
        raw_path.name,
    )

    print(
        "CTD report metadata:",
        report_created,
    )

    print()

    print(
        "Human curated rows:",
        f"{human_rows:,}",
    )

    print(
        "Mapped human rows:",
        f"{mapped_human_rows:,}",
    )

    print(
        "Rows with gene outside "
        "Core V1:",
        f"{unmapped_gene_rows:,}",
    )

    print(
        "Rows rejected as complex:",
        f"{complex_rows_rejected:,}",
    )

    print(
        "Rows rejected as cotreatment:",
        f"{cotreatment_rows_rejected:,}",
    )

    print(
        "Rows rejected as no-effect:",
        f"{no_effect_rows_rejected:,}",
    )

    print()

    print(
        "Strict retained curated rows:",
        f"{strict_curated_rows:,}",
    )

    print(
        "Chemical nodes:",
        f"{len(chemical_rows):,}",
    )

    print(
        "Chemical nodes without CAS RN:",
        f"{missing_cas_chemicals:,}",
    )

    print(
        "Chemical-Gene edges:",
        f"{len(edge_rows):,}",
    )

    print(
        "Core genes with >=1 "
        "Chemical-Gene edge:",
        f"{len(active_core_genes):,}",
    )

    print(
        "Core gene coverage:",
        f"{coverage:.2f}%",
    )

    print()

    print(
        "Chemicals written:",
        CHEMICALS_OUTPUT.relative_to(
            PROJECT_ROOT
        ),
    )

    print(
        "Chemical-Gene layer written:",
        EDGES_OUTPUT.relative_to(
            PROJECT_ROOT
        ),
    )


if __name__ == "__main__":
    main()