from __future__ import annotations

import csv
import gzip
import hashlib
import re
import statistics
import tomllib
from collections import Counter
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

TISSUES_OUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "nodes"
    / "tissues.tsv"
)

EDGES_OUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "edges"
    / "gene_tissue.tsv"
)


HGNC_RE = re.compile(
    r"^HGNC:(\d+)$"
)

UBERON_RE = re.compile(
    r"^UBERON:(\d+)$"
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
        "Verified SHA-256:",
        path.name,
    )


def hgnc_sort_key(
    value: str,
) -> int:
    match = HGNC_RE.fullmatch(
        value
    )

    if not match:
        raise ValueError(
            f"Invalid HGNC ID: {value}"
        )

    return int(
        match.group(1)
    )


def uberon_sort_key(
    value: str,
) -> int:
    match = UBERON_RE.fullmatch(
        value
    )

    if not match:
        raise ValueError(
            f"Invalid UBERON ID: {value}"
        )

    return int(
        match.group(1)
    )


def load_core_genes(
    path: Path,
) -> tuple[
    set[str],
    dict[str, str],
]:
    core_ids: set[str] = set()

    ensembl_to_hgnc: dict[
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
            "ensembl_gene_id",
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

            ensembl_id = row[
                "ensembl_gene_id"
            ].strip()

            if (
                not ensembl_id
                or ensembl_id == "NA"
            ):
                continue

            previous = (
                ensembl_to_hgnc.get(
                    ensembl_id
                )
            )

            if (
                previous is not None
                and previous != hgnc_id
            ):
                raise ValueError(
                    "Ensembl Gene ID maps "
                    "to multiple Core genes: "
                    f"{ensembl_id}"
                )

            ensembl_to_hgnc[
                ensembl_id
            ] = hgnc_id

    return (
        core_ids,
        ensembl_to_hgnc,
    )


def write_tissues(
    output_path: Path,
    tissue_names: dict[
        str,
        str,
    ],
    retained_tissues: set[str],
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
            ]
        )

        for tissue_id in sorted(
            retained_tissues,
            key=uberon_sort_key,
        ):
            writer.writerow(
                [
                    tissue_id,
                    tissue_names[
                        tissue_id
                    ],
                ]
            )


def write_edges(
    output_path: Path,
    selected_rows: dict,
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
            uberon_sort_key(
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
                "bgee_gene_id",
                "call_quality",
                "fdr",
                "expression_score",
                "expression_rank",
                "self_observation_count",
                "descendant_observation_count",
            ]
        )

        for pair in sorted(
            selected_rows,
            key=edge_sort_key,
        ):
            gene_id, tissue_id = pair

            item = selected_rows[
                pair
            ]

            writer.writerow(
                [
                    gene_id,
                    tissue_id,
                    relation,
                    item[
                        "bgee_gene_id"
                    ],
                    item[
                        "call_quality"
                    ],
                    item[
                        "fdr_text"
                    ],
                    item[
                        "score_text"
                    ],
                    item[
                        "rank_text"
                    ],
                    item[
                        "self_observation_count"
                    ],
                    item[
                        "descendant_observation_count"
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
        "gene_tissue_network"
    ]

    expected = policy[
        "expected_build"
    ]

    source = sources[
        "sources"
    ][
        "bgee"
    ]

    # ----------------------------------------------
    # Configuration invariants
    # ----------------------------------------------

    if policy["source"] != "Bgee":
        raise ValueError(
            "gene_tissue_network.source "
            "must be Bgee"
        )

    if (
        str(
            policy[
                "source_version"
            ]
        )
        != str(
            source[
                "version"
            ]
        )
    ):
        raise ValueError(
            "Bgee version mismatch "
            "between dataset and "
            "source configs"
        )

    if (
        policy[
            "anatomical_namespace"
        ]
        != "UBERON"
    ):
        raise ValueError(
            "Tissue V1 requires "
            "UBERON-only anatomy"
        )

    if (
        policy[
            "selection_strategy"
        ]
        !=
        "top_k_per_gene_by_expression_rank"
    ):
        raise ValueError(
            "Unexpected Bgee "
            "selection strategy"
        )

    if not policy[
        "lower_expression_rank_is_better"
    ]:
        raise ValueError(
            "Bgee expression rank "
            "must be interpreted as "
            "lower-is-better"
        )

    if (
        policy[
            "tie_policy"
        ]
        != "retain_all_boundary_ties"
    ):
        raise ValueError(
            "Unexpected Top-K tie policy"
        )

    if (
        source[
            "condition_scope"
        ]
        != "anatomical_entities_only"
    ):
        raise ValueError(
            "Bgee source must use "
            "anatomical-entities-only calls"
        )

    if (
        source[
            "file_format"
        ]
        != "advanced"
    ):
        raise ValueError(
            "Bgee source must use "
            "advanced calls file"
        )

    # ----------------------------------------------
    # Raw source
    # ----------------------------------------------

    raw_path = (
        PROJECT_ROOT
        / source[
            "local_path"
        ]
    )

    verify_sha256(
        raw_path,
        source[
            "sha256"
        ],
    )

    # ----------------------------------------------
    # Frozen Core genes
    # ----------------------------------------------

    (
        core_ids,
        ensembl_to_hgnc,
    ) = load_core_genes(
        GENES_FILE
    )

    # ----------------------------------------------
    # Policy
    # ----------------------------------------------

    required_expression = policy[
        "expression_status"
    ]

    required_quality = policy[
        "call_quality"
    ]

    top_k = int(
        policy[
            "top_k"
        ]
    )

    if top_k < 1:
        raise ValueError(
            "top_k must be >= 1"
        )

    excluded_tissues = set(
        policy[
            "exclude_anatomical_ids"
        ]
    )

    for tissue_id in excluded_tissues:
        uberon_sort_key(
            tissue_id
        )

    maximum_allowed_fdr = float(
        expected[
            "maximum_fdr"
        ]
    )

    # ----------------------------------------------
    # Output accumulators
    #
    # Only final Top-K rows are kept in memory.
    # ----------------------------------------------

    selected_rows: dict[
        tuple[str, str],
        dict,
    ] = {}

    tissue_names: dict[
        str,
        str,
    ] = {}

    seen_source_gene_ids: set[
        str
    ] = set()

    current_bgee_gene_id = None
    current_hgnc_id = None
    current_rows: list[dict] = []

    candidate_rows = 0
    excluded_tissue_rows = 0
    mapped_source_genes = set()

    genes_below_k = 0
    genes_exact_k = 0
    genes_above_k_due_to_ties = 0

    # ----------------------------------------------
    # Per-gene Top-K selector
    # ----------------------------------------------

    def process_current_gene() -> None:
        nonlocal genes_below_k
        nonlocal genes_exact_k
        nonlocal genes_above_k_due_to_ties

        if (
            current_bgee_gene_id is None
            or current_hgnc_id is None
            or not current_rows
        ):
            return

        # One call per gene × anatomical entity
        # is expected in an anatomy-only Bgee file.
        tissue_ids = [
            item["tissue_id"]
            for item in current_rows
        ]

        if (
            len(tissue_ids)
            != len(set(tissue_ids))
        ):
            raise ValueError(
                "Duplicate Bgee Gene-Tissue "
                "candidate rows for "
                f"{current_bgee_gene_id}"
            )

        ordered = sorted(
            current_rows,
            key=lambda item: (
                item[
                    "rank"
                ],
                item[
                    "fdr"
                ],
                item[
                    "tissue_id"
                ],
            ),
        )

        if len(ordered) <= top_k:
            retained = ordered
        else:
            boundary_rank = ordered[
                top_k - 1
            ][
                "rank"
            ]

            retained = [
                item
                for item
                in ordered
                if (
                    item[
                        "rank"
                    ]
                    <= boundary_rank
                )
            ]

        retained_count = len(
            retained
        )

        if retained_count < top_k:
            genes_below_k += 1

        elif retained_count == top_k:
            genes_exact_k += 1

        else:
            genes_above_k_due_to_ties += 1

        for item in retained:

            pair = (
                current_hgnc_id,
                item[
                    "tissue_id"
                ],
            )

            if pair in selected_rows:
                raise ValueError(
                    "Duplicate final "
                    "Gene-Tissue pair: "
                    f"{pair}"
                )

            selected_rows[
                pair
            ] = item

    # ----------------------------------------------
    # Stream Bgee source
    # ----------------------------------------------

    with gzip.open(
        raw_path,
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
            quotechar='"',
        )

        required_columns = {
            "Gene ID",
            "Anatomical entity ID",
            "Anatomical entity name",
            "Expression",
            "Call quality",
            "FDR",
            "Expression score",
            "Expression rank",
            "Including observed data",
            "Self observation count",
            "Descendant observation count",
        }

        if reader.fieldnames is None:
            raise ValueError(
                "Bgee file has no header"
            )

        missing = (
            required_columns
            - set(
                reader.fieldnames
            )
        )

        if missing:
            raise ValueError(
                "Bgee file missing "
                "required columns: "
                f"{sorted(missing)}"
            )

        print(
            "Bgee columns:",
            len(
                reader.fieldnames
            ),
        )

        for line_number, row in enumerate(
            reader,
            start=2,
        ):

            bgee_gene_id = row[
                "Gene ID"
            ].strip()

            if not bgee_gene_id:
                raise ValueError(
                    "Missing Bgee Gene ID "
                    f"at line {line_number}"
                )

            # --------------------------------------
            # Source file grouping guard
            # --------------------------------------

            if (
                current_bgee_gene_id
                is None
            ):
                current_bgee_gene_id = (
                    bgee_gene_id
                )

                base_ensembl = (
                    bgee_gene_id.split(
                        ".",
                        1,
                    )[0]
                )

                current_hgnc_id = (
                    ensembl_to_hgnc.get(
                        base_ensembl
                    )
                )

            elif (
                bgee_gene_id
                != current_bgee_gene_id
            ):
                process_current_gene()

                seen_source_gene_ids.add(
                    current_bgee_gene_id
                )

                if (
                    bgee_gene_id
                    in seen_source_gene_ids
                ):
                    raise ValueError(
                        "Bgee source is not "
                        "grouped by Gene ID: "
                        f"{bgee_gene_id}"
                    )

                current_rows = []

                current_bgee_gene_id = (
                    bgee_gene_id
                )

                base_ensembl = (
                    bgee_gene_id.split(
                        ".",
                        1,
                    )[0]
                )

                current_hgnc_id = (
                    ensembl_to_hgnc.get(
                        base_ensembl
                    )
                )

            # --------------------------------------
            # Only Core genes
            # --------------------------------------

            if current_hgnc_id is None:
                continue

            mapped_source_genes.add(
                current_hgnc_id
            )

            # --------------------------------------
            # V1 expression policy
            # --------------------------------------

            if (
                row[
                    "Expression"
                ].strip()
                != required_expression
            ):
                continue

            if (
                row[
                    "Call quality"
                ].strip()
                != required_quality
            ):
                continue

            tissue_id = row[
                "Anatomical entity ID"
            ].strip()

            # Pure UBERON only.
            if not UBERON_RE.fullmatch(
                tissue_id
            ):
                continue

            if tissue_id in excluded_tissues:
                excluded_tissue_rows += 1
                continue

            if (
                row[
                    "Including observed data"
                ].strip()
                != "yes"
            ):
                raise ValueError(
                    "Selected Bgee call "
                    "without observed data: "
                    f"{bgee_gene_id} / "
                    f"{tissue_id}"
                )

            tissue_name = row[
                "Anatomical entity name"
            ].strip()

            if not tissue_name:
                raise ValueError(
                    "Missing tissue name "
                    f"for {tissue_id}"
                )

            previous_name = (
                tissue_names.get(
                    tissue_id
                )
            )

            if (
                previous_name is not None
                and previous_name
                != tissue_name
            ):
                raise ValueError(
                    "Tissue has multiple names: "
                    f"{tissue_id}: "
                    f"{previous_name!r} vs "
                    f"{tissue_name!r}"
                )

            tissue_names[
                tissue_id
            ] = tissue_name

            fdr_text = row[
                "FDR"
            ].strip()

            score_text = row[
                "Expression score"
            ].strip()

            rank_text = row[
                "Expression rank"
            ].strip()

            try:
                fdr = float(
                    fdr_text
                )

                score = float(
                    score_text
                )

                rank = float(
                    rank_text
                )

                self_count = int(
                    row[
                        "Self observation count"
                    ]
                )

                descendant_count = int(
                    row[
                        "Descendant observation count"
                    ]
                )

            except ValueError as error:
                raise ValueError(
                    "Invalid numeric Bgee "
                    "field at line "
                    f"{line_number}"
                ) from error

            if fdr > maximum_allowed_fdr:
                raise ValueError(
                    "Gold call exceeds "
                    "configured maximum FDR: "
                    f"{bgee_gene_id} / "
                    f"{tissue_id} / "
                    f"{fdr}"
                )

            if self_count < 1:
                raise ValueError(
                    "Selected Bgee call has "
                    "no self observation: "
                    f"{bgee_gene_id} / "
                    f"{tissue_id}"
                )

            if descendant_count < 0:
                raise ValueError(
                    "Negative descendant "
                    "observation count"
                )

            candidate_rows += 1

            current_rows.append(
                {
                    "bgee_gene_id":
                        bgee_gene_id,
                    "tissue_id":
                        tissue_id,
                    "call_quality":
                        required_quality,
                    "fdr":
                        fdr,
                    "fdr_text":
                        fdr_text,
                    "score":
                        score,
                    "score_text":
                        score_text,
                    "rank":
                        rank,
                    "rank_text":
                        rank_text,
                    "self_observation_count":
                        self_count,
                    "descendant_observation_count":
                        descendant_count,
                }
            )

    # Flush final source gene.
    process_current_gene()

    # ----------------------------------------------
    # Final topology
    # ----------------------------------------------

    active_genes = {
        gene_id
        for (
            gene_id,
            _tissue_id,
        )
        in selected_rows
    }

    retained_tissues = {
        tissue_id
        for (
            _gene_id,
            tissue_id,
        )
        in selected_rows
    }

    gene_degrees = Counter(
        gene_id
        for (
            gene_id,
            _tissue_id,
        )
        in selected_rows
    )

    tissue_degrees = Counter(
        tissue_id
        for (
            _gene_id,
            tissue_id,
        )
        in selected_rows
    )

    # ----------------------------------------------
    # Regression guards
    # ----------------------------------------------

    if (
        len(retained_tissues)
        != int(
            expected[
                "tissues"
            ]
        )
    ):
        raise ValueError(
            "Tissue regression failure: "
            f"expected "
            f"{expected['tissues']}, "
            f"observed "
            f"{len(retained_tissues)}"
        )

    if (
        len(selected_rows)
        != int(
            expected[
                "edges"
            ]
        )
    ):
        raise ValueError(
            "Gene-Tissue edge "
            "regression failure: "
            f"expected "
            f"{expected['edges']}, "
            f"observed "
            f"{len(selected_rows)}"
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
        min(
            gene_degrees.values()
        )
        != int(
            expected[
                "minimum_gene_degree"
            ]
        )
    ):
        raise ValueError(
            "Minimum gene-degree "
            "regression failure"
        )

    if (
        statistics.median(
            gene_degrees.values()
        )
        != float(
            expected[
                "median_gene_degree"
            ]
        )
    ):
        raise ValueError(
            "Median gene-degree "
            "regression failure"
        )

    if (
        max(
            gene_degrees.values()
        )
        != int(
            expected[
                "maximum_gene_degree"
            ]
        )
    ):
        raise ValueError(
            "Maximum gene-degree "
            "regression failure"
        )

    if (
        max(
            tissue_degrees.values()
        )
        != int(
            expected[
                "maximum_tissue_degree"
            ]
        )
    ):
        raise ValueError(
            "Maximum tissue-degree "
            "regression failure"
        )

    observed_max_fdr = max(
        item[
            "fdr"
        ]
        for item
        in selected_rows.values()
    )

    if (
        observed_max_fdr
        > maximum_allowed_fdr
    ):
        raise ValueError(
            "Maximum FDR regression "
            "failure"
        )

    observed_coverage = (
        100.0
        * len(active_genes)
        / len(core_ids)
    )

    expected_coverage = float(
        expected[
            "core_gene_coverage_percent"
        ]
    )

    if (
        abs(
            observed_coverage
            - expected_coverage
        )
        > 0.01
    ):
        raise ValueError(
            "Core coverage regression "
            "failure: expected "
            f"{expected_coverage:.2f}%, "
            f"observed "
            f"{observed_coverage:.2f}%"
        )

    # ----------------------------------------------
    # Final invariants
    # ----------------------------------------------

    if not active_genes <= core_ids:
        raise ValueError(
            "Gene-Tissue edges contain "
            "genes outside Core"
        )

    if (
        retained_tissues
        & excluded_tissues
    ):
        raise ValueError(
            "Excluded anatomy appears "
            "in final Tissue layer"
        )

    if any(
        not UBERON_RE.fullmatch(
            tissue_id
        )
        for tissue_id
        in retained_tissues
    ):
        raise ValueError(
            "Non-UBERON node appears "
            "in Tissue layer"
        )

    if any(
        degree > (
            top_k
            + 20
        )
        for degree
        in gene_degrees.values()
    ):
        raise ValueError(
            "Unexpectedly large "
            "Top-K tie expansion"
        )

    # ----------------------------------------------
    # Write outputs
    # ----------------------------------------------

    relation = dataset[
        "edge_types"
    ][
        "gene_tissue"
    ][
        "relation"
    ]

    write_tissues(
        TISSUES_OUT,
        tissue_names,
        retained_tissues,
    )

    write_edges(
        EDGES_OUT,
        selected_rows,
        relation,
    )

    # ----------------------------------------------
    # Report
    # ----------------------------------------------

    print()

    print(
        "Bgee Tissue V1 build complete"
    )

    print(
        "Mapped Bgee Core genes:",
        f"{len(mapped_source_genes):,}",
    )

    print(
        "Gold-present UBERON "
        "candidate rows:",
        f"{candidate_rows:,}",
    )

    print(
        "Excluded generic-tissue rows:",
        f"{excluded_tissue_rows:,}",
    )

    print(
        "Retained tissue nodes:",
        f"{len(retained_tissues):,}",
    )

    print(
        "Gene-Tissue edges:",
        f"{len(selected_rows):,}",
    )

    print(
        "Active Core genes:",
        (
            f"{len(active_genes):,} / "
            f"{len(core_ids):,} "
            f"({observed_coverage:.2f}%)"
        ),
    )

    print(
        "Genes with < Top-K tissues:",
        f"{genes_below_k:,}",
    )

    print(
        "Genes with exactly Top-K:",
        f"{genes_exact_k:,}",
    )

    print(
        "Genes exceeding Top-K "
        "because of rank ties:",
        f"{genes_above_k_due_to_ties:,}",
    )

    print(
        "Gene degree:",
        (
            f"min={min(gene_degrees.values())}, "
            f"median="
            f"{statistics.median(gene_degrees.values())}, "
            f"max={max(gene_degrees.values())}"
        ),
    )

    print(
        "Maximum tissue degree:",
        f"{max(tissue_degrees.values()):,}",
    )

    print(
        "Maximum retained FDR:",
        observed_max_fdr,
    )

    print(
        "Wrote:",
        TISSUES_OUT.relative_to(
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