from __future__ import annotations

import csv
import hashlib
import re
import tomllib
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATASET_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "dataset_v1.toml"
)

SOURCES_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "sources_v1.toml"
)

GENES_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "nodes"
    / "genes.tsv"
)

PATHWAYS_OUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "nodes"
    / "pathways.tsv"
)

EDGES_OUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "edges"
    / "gene_pathway.tsv"
)


REACTOME_ID_RE = re.compile(
    r"^R-HSA-\d+$"
)

HGNC_ID_RE = re.compile(
    r"^HGNC:(\d+)$"
)


def load_toml(
    path: Path,
) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def sha256_file(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def verify_sha256(
    path: Path,
    expected: str,
) -> None:
    observed = sha256_file(
        path
    )

    if (
        observed.lower()
        != expected.lower()
    ):
        raise ValueError(
            "SHA-256 mismatch for "
            f"{path.name}: "
            f"expected {expected}, "
            f"observed {observed}"
        )

    print(
        f"Verified SHA-256: "
        f"{path.name}"
    )


def hgnc_sort_key(
    value: str,
) -> int:
    match = HGNC_ID_RE.fullmatch(
        value
    )

    if not match:
        raise ValueError(
            f"Invalid HGNC ID: {value}"
        )

    return int(
        match.group(1)
    )


def reactome_sort_key(
    value: str,
) -> int:
    if not REACTOME_ID_RE.fullmatch(
        value
    ):
        raise ValueError(
            "Invalid human Reactome ID: "
            f"{value}"
        )

    return int(
        value.rsplit(
            "-",
            1,
        )[1]
    )


def load_core_genes(
    path: Path,
) -> tuple[
    set[str],
    dict[str, str],
]:
    core_ids: set[str] = set()

    ncbi_to_hgnc: dict[
        str,
        str,
    ] = {}

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "canonical_id",
            "ncbi_gene_id",
        }

        if (
            reader.fieldnames is None
            or not required.issubset(
                reader.fieldnames
            )
        ):
            raise ValueError(
                "genes.tsv missing "
                "required columns: "
                f"{sorted(required)}"
            )

        for row in reader:

            hgnc_id = row[
                "canonical_id"
            ].strip()

            hgnc_sort_key(
                hgnc_id
            )

            if hgnc_id in core_ids:
                raise ValueError(
                    "Duplicate Core HGNC ID: "
                    f"{hgnc_id}"
                )

            core_ids.add(
                hgnc_id
            )

            ncbi_id = row[
                "ncbi_gene_id"
            ].strip()

            if (
                not ncbi_id
                or ncbi_id == "NA"
            ):
                continue

            previous = (
                ncbi_to_hgnc.get(
                    ncbi_id
                )
            )

            if (
                previous is not None
                and previous != hgnc_id
            ):
                raise ValueError(
                    "NCBI Gene ID maps "
                    "to multiple Core genes: "
                    f"{ncbi_id} -> "
                    f"{previous}, {hgnc_id}"
                )

            ncbi_to_hgnc[
                ncbi_id
            ] = hgnc_id

    return (
        core_ids,
        ncbi_to_hgnc,
    )


def load_human_pathways(
    path: Path,
) -> dict[str, str]:
    pathways: dict[
        str,
        str,
    ] = {}

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:

        for (
            line_number,
            line,
        ) in enumerate(
            handle,
            start=1,
        ):

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            parts = line.split(
                "\t"
            )

            if len(parts) != 3:
                raise ValueError(
                    "Unexpected "
                    "ReactomePathways row "
                    f"at line {line_number}: "
                    f"{len(parts)} columns"
                )

            (
                pathway_id,
                name,
                species,
            ) = parts

            if (
                species
                != "Homo sapiens"
            ):
                continue

            reactome_sort_key(
                pathway_id
            )

            if not name:
                raise ValueError(
                    "Missing pathway name "
                    f"for {pathway_id}"
                )

            if pathway_id in pathways:
                raise ValueError(
                    "Duplicate human "
                    "pathway ID: "
                    f"{pathway_id}"
                )

            pathways[
                pathway_id
            ] = name

    return pathways


def load_disease_pathways(
    path: Path,
    human_pathways: dict[
        str,
        str,
    ],
) -> set[str]:
    disease_ids: set[str] = set()

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:

        for (
            line_number,
            line,
        ) in enumerate(
            handle,
            start=1,
        ):

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            parts = line.split(
                "\t"
            )

            if len(parts) != 2:
                raise ValueError(
                    "Unexpected "
                    "HumanDiseasePathways row "
                    f"at line {line_number}: "
                    f"{len(parts)} columns"
                )

            pathway_id, _ = parts

            if (
                pathway_id
                not in human_pathways
            ):
                raise ValueError(
                    "Disease pathway absent "
                    "from human master: "
                    f"{pathway_id}"
                )

            disease_ids.add(
                pathway_id
            )

    return disease_ids


def build_pairs(
    mapping_path: Path,
    ncbi_to_hgnc: dict[
        str,
        str,
    ],
    human_pathways: dict[
        str,
        str,
    ],
    allowed_evidence: list[str],
):
    allowed = set(
        allowed_evidence
    )

    pair_evidence = defaultdict(
        set
    )

    pair_support = Counter()

    human_rows = 0
    mapped_rows = 0

    unmapped_ncbi_ids: set[
        str
    ] = set()

    with mapping_path.open(
        encoding="utf-8",
        newline="",
    ) as handle:

        for (
            line_number,
            line,
        ) in enumerate(
            handle,
            start=1,
        ):

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            parts = line.split(
                "\t"
            )

            if len(parts) != 6:
                raise ValueError(
                    "Unexpected "
                    "NCBI2Reactome row "
                    f"at line {line_number}: "
                    f"{len(parts)} columns"
                )

            (
                ncbi_id,
                pathway_id,
                _url,
                _mapping_name,
                evidence,
                species,
            ) = parts

            if (
                species
                != "Homo sapiens"
            ):
                continue

            human_rows += 1

            if (
                evidence
                not in allowed
            ):
                continue

            hgnc_id = (
                ncbi_to_hgnc.get(
                    ncbi_id
                )
            )

            if hgnc_id is None:
                unmapped_ncbi_ids.add(
                    ncbi_id
                )
                continue

            if (
                pathway_id
                not in human_pathways
            ):
                raise ValueError(
                    "Human mapping references "
                    "pathway absent from "
                    "human master: "
                    f"{pathway_id}"
                )

            mapped_rows += 1

            pair = (
                hgnc_id,
                pathway_id,
            )

            pair_evidence[
                pair
            ].add(
                evidence
            )

            pair_support[
                pair
            ] += 1

    return (
        pair_evidence,
        pair_support,
        human_rows,
        mapped_rows,
        unmapped_ncbi_ids,
    )


def write_pathways(
    output_path: Path,
    retained_pathways: set[str],
    pathway_names: dict[
        str,
        str,
    ],
    disease_pathways: set[str],
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "canonical_id",
                "name",
                "is_disease_pathway",
            ]
        )

        for pathway_id in sorted(
            retained_pathways,
            key=reactome_sort_key,
        ):
            writer.writerow(
                [
                    pathway_id,
                    pathway_names[
                        pathway_id
                    ],
                    (
                        "true"
                        if pathway_id
                        in disease_pathways
                        else "false"
                    ),
                ]
            )


def write_edges(
    output_path: Path,
    retained_pairs,
    pair_evidence,
    pair_support,
    evidence_order: list[str],
    relation: str,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    def edge_sort_key(
        pair,
    ):
        return (
            hgnc_sort_key(
                pair[0]
            ),
            reactome_sort_key(
                pair[1]
            ),
        )

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "source_id",
                "target_id",
                "relation",
                "evidence_codes",
                "supporting_mapping_rows",
            ]
        )

        for pair in sorted(
            retained_pairs,
            key=edge_sort_key,
        ):
            (
                gene_id,
                pathway_id,
            ) = pair

            evidence_codes = "|".join(
                code
                for code
                in evidence_order
                if code
                in pair_evidence[
                    pair
                ]
            )

            writer.writerow(
                [
                    gene_id,
                    pathway_id,
                    relation,
                    evidence_codes,
                    pair_support[
                        pair
                    ],
                ]
            )


def main() -> None:
    dataset = load_toml(
        DATASET_CONFIG
    )

    sources = load_toml(
        SOURCES_CONFIG
    )

    policy = dataset[
        "gene_pathway_network"
    ]

    expected = policy[
        "expected_build"
    ]

    source = sources[
        "sources"
    ][
        "reactome"
    ]

    # ----------------------------------------------
    # Configuration invariants
    # ----------------------------------------------

    if (
        policy["source"]
        != "Reactome"
    ):
        raise ValueError(
            "gene_pathway_network.source "
            "must be Reactome"
        )

    if (
        str(
            policy[
                "source_version"
            ]
        )
        != str(
            source["version"]
        )
    ):
        raise ValueError(
            "Reactome version mismatch "
            "between dataset and "
            "source configs"
        )

    if (
        policy[
            "primary_mapping"
        ]
        != "NCBI2Reactome"
    ):
        raise ValueError(
            "Reactome V1 primary "
            "mapping must be "
            "NCBI2Reactome"
        )

    if (
        policy[
            "pathway_level"
        ]
        != "lowest_level"
    ):
        raise ValueError(
            "Reactome V1 must use "
            "lowest-level pathway "
            "mappings"
        )

    if not policy[
        "include_disease_pathways"
    ]:
        raise ValueError(
            "Reactome V1 is configured "
            "to retain disease pathways"
        )

    # ----------------------------------------------
    # Raw source paths
    # ----------------------------------------------

    ncbi_map = (
        PROJECT_ROOT
        / source[
            "ncbi_mapping_local_path"
        ]
    )

    ensembl_map = (
        PROJECT_ROOT
        / source[
            "ensembl_mapping_local_path"
        ]
    )

    pathways_file = (
        PROJECT_ROOT
        / source[
            "pathways_local_path"
        ]
    )

    disease_file = (
        PROJECT_ROOT
        / source[
            "disease_pathways_local_path"
        ]
    )

    # ----------------------------------------------
    # Immutable source verification
    # ----------------------------------------------

    verify_sha256(
        ncbi_map,
        source[
            "ncbi_mapping_sha256"
        ],
    )

    verify_sha256(
        ensembl_map,
        source[
            "ensembl_mapping_sha256"
        ],
    )

    verify_sha256(
        pathways_file,
        source[
            "pathways_sha256"
        ],
    )

    verify_sha256(
        disease_file,
        source[
            "disease_pathways_sha256"
        ],
    )

    # ----------------------------------------------
    # Frozen gene universe
    # ----------------------------------------------

    (
        core_ids,
        ncbi_to_hgnc,
    ) = load_core_genes(
        GENES_FILE
    )

    # ----------------------------------------------
    # Reactome pathway metadata
    # ----------------------------------------------

    pathway_names = (
        load_human_pathways(
            pathways_file
        )
    )

    disease_pathways = (
        load_disease_pathways(
            disease_file,
            pathway_names,
        )
    )

    # ----------------------------------------------
    # Build canonical Gene-Pathway pairs
    # ----------------------------------------------

    evidence_order = list(
        policy[
            "include_evidence_codes"
        ]
    )

    if not evidence_order:
        raise ValueError(
            "include_evidence_codes "
            "must not be empty"
        )

    (
        pair_evidence,
        pair_support,
        human_rows,
        mapped_rows,
        unmapped_ncbi_ids,
    ) = build_pairs(
        ncbi_map,
        ncbi_to_hgnc,
        pathway_names,
        evidence_order,
    )

    # ----------------------------------------------
    # Pathway-size filtering
    # ----------------------------------------------

    pathway_sizes = Counter(
        pathway_id
        for _,
        pathway_id
        in pair_evidence
    )

    minimum_size = int(
        policy[
            "minimum_core_gene_count"
        ]
    )

    maximum_size = int(
        policy[
            "maximum_core_gene_count"
        ]
    )

    retained_pathways = {
        pathway_id
        for (
            pathway_id,
            size,
        ) in pathway_sizes.items()
        if (
            size >= minimum_size
            and (
                maximum_size == 0
                or size <= maximum_size
            )
        )
    }

    retained_pairs = {
        pair
        for pair
        in pair_evidence
        if pair[1]
        in retained_pathways
    }

    active_genes = {
        gene_id
        for (
            gene_id,
            _,
        ) in retained_pairs
    }

    retained_disease = (
        retained_pathways
        & disease_pathways
    )

    retained_sizes = Counter(
        pathway_id
        for _,
        pathway_id
        in retained_pairs
    )

    max_pathway_size = max(
        retained_sizes.values(),
        default=0,
    )

    # ----------------------------------------------
    # Regression guards
    # ----------------------------------------------

    if (
        len(retained_pathways)
        != int(
            expected[
                "pathways"
            ]
        )
    ):
        raise ValueError(
            "Pathway regression failure: "
            f"expected "
            f"{expected['pathways']}, "
            f"observed "
            f"{len(retained_pathways)}"
        )

    if (
        len(retained_disease)
        != int(
            expected[
                "disease_pathways"
            ]
        )
    ):
        raise ValueError(
            "Disease-pathway regression "
            "failure: expected "
            f"{expected['disease_pathways']}, "
            f"observed "
            f"{len(retained_disease)}"
        )

    if (
        len(retained_pairs)
        != int(
            expected[
                "edges"
            ]
        )
    ):
        raise ValueError(
            "Edge regression failure: "
            f"expected "
            f"{expected['edges']}, "
            f"observed "
            f"{len(retained_pairs)}"
        )

    if (
        len(active_genes)
        != int(
            expected[
                "active_core_genes"
            ]
        )
    ):
        raise ValueError(
            "Active-gene regression "
            "failure: expected "
            f"{expected['active_core_genes']}, "
            f"observed "
            f"{len(active_genes)}"
        )

    if (
        max_pathway_size
        != int(
            expected[
                "maximum_pathway_size"
            ]
        )
    ):
        raise ValueError(
            "Maximum pathway-size "
            "regression failure: "
            f"expected "
            f"{expected['maximum_pathway_size']}, "
            f"observed "
            f"{max_pathway_size}"
        )

    # ----------------------------------------------
    # Final invariants
    # ----------------------------------------------

    if not active_genes <= core_ids:
        raise ValueError(
            "Gene-Pathway edges "
            "contain genes outside Core"
        )

    if not (
        retained_pathways
        <= set(pathway_names)
    ):
        raise ValueError(
            "Gene-Pathway edges "
            "contain unknown pathways"
        )

    if any(
        pair_support[pair] < 1
        for pair
        in retained_pairs
    ):
        raise ValueError(
            "Invalid "
            "supporting_mapping_rows "
            "value"
        )

    if any(
        not pair_evidence[pair]
        for pair
        in retained_pairs
    ):
        raise ValueError(
            "Gene-Pathway edge "
            "without evidence code"
        )

    # ----------------------------------------------
    # Write outputs
    # ----------------------------------------------

    relation = dataset[
        "edge_types"
    ][
        "gene_pathway"
    ][
        "relation"
    ]

    write_pathways(
        PATHWAYS_OUT,
        retained_pathways,
        pathway_names,
        disease_pathways,
    )

    write_edges(
        EDGES_OUT,
        retained_pairs,
        pair_evidence,
        pair_support,
        evidence_order,
        relation,
    )

    # ----------------------------------------------
    # Build report
    # ----------------------------------------------

    coverage = (
        100.0
        * len(active_genes)
        / len(core_ids)
    )

    print()

    print(
        "Reactome V1 build complete"
    )

    print(
        "Human NCBI2Reactome rows:",
        f"{human_rows:,}",
    )

    print(
        "Mapped human rows:",
        f"{mapped_rows:,}",
    )

    print(
        "Unmapped human NCBI IDs:",
        f"{len(unmapped_ncbi_ids):,}",
    )

    print(
        "Unique mapped Gene-Pathway "
        "pairs before size filter:",
        f"{len(pair_evidence):,}",
    )

    print(
        "Retained pathways:",
        f"{len(retained_pathways):,}",
    )

    print(
        "Retained disease pathways:",
        f"{len(retained_disease):,}",
    )

    print(
        "Gene-Pathway edges:",
        f"{len(retained_pairs):,}",
    )

    print(
        "Active Core genes:",
        (
            f"{len(active_genes):,} / "
            f"{len(core_ids):,} "
            f"({coverage:.2f}%)"
        ),
    )

    print(
        "Maximum retained "
        "pathway size:",
        f"{max_pathway_size:,}",
    )

    print(
        "Wrote:",
        PATHWAYS_OUT.relative_to(
            PROJECT_ROOT
        ),
    )

    print(
        "Wrote:",
        EDGES_OUT.relative_to(
            PROJECT_ROOT
        ),
    )


if __name__ == "__main__":
    main()