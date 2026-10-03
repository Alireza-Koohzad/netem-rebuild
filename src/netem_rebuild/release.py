from __future__ import annotations

import csv
import hashlib
import json
import tomllib
from datetime import date, datetime
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_ROOT = (
    PROJECT_ROOT
    / "data"
)

OUTPUT = (
    DATA_ROOT
    / "datapackage.json"
)

DATASET_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "dataset_v1.toml"
)

CONTRACT_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "data_contract_v1.toml"
)

SOURCES_CONFIG = (
    PROJECT_ROOT
    / "configs"
    / "sources_v1.toml"
)


EXPECTED_TOTAL_NODES = 30_832
EXPECTED_TOTAL_EDGES = 1_614_468


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


def json_safe(
    value: Any,
) -> Any:
    if isinstance(
        value,
        (date, datetime),
    ):
        return value.isoformat()

    if isinstance(
        value,
        Path,
    ):
        return value.as_posix()

    if isinstance(
        value,
        dict,
    ):
        return {
            str(key): json_safe(item)
            for key, item
            in value.items()
        }

    if isinstance(
        value,
        (list, tuple),
    ):
        return [
            json_safe(item)
            for item in value
        ]

    return value


def inspect_tsv(
    path: Path,
) -> tuple[
    list[str],
    int,
]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.reader(
            handle,
            delimiter="\t",
        )

        try:
            header = next(
                reader
            )
        except StopIteration as error:
            raise ValueError(
                f"Empty TSV file: {path}"
            ) from error

        row_count = sum(
            1
            for _ in reader
        )

    return (
        header,
        row_count,
    )


def field(
    name: str,
    field_type: str = "string",
) -> dict:
    return {
        "name": name,
        "type": field_type,
    }


RESOURCE_SPECS = [
    {
        "name": "genes",
        "path": "processed/nodes/genes.tsv",
        "kind": "node",
        "node_type": "gene",
        "expected_rows": 19_297,
        "fields": [
            field("canonical_id"),
            field("symbol"),
            field("name"),
            field("ensembl_gene_id"),
            field("ncbi_gene_id"),
        ],
        "primary_key": "canonical_id",
        "foreign_keys": [],
    },
    {
        "name": "chemicals",
        "path": "processed/nodes/chemicals.tsv",
        "kind": "node",
        "node_type": "chemical",
        "expected_rows": 9_601,
        "fields": [
            field("canonical_id"),
            field("name"),
            field("cas_rn"),
        ],
        "primary_key": "canonical_id",
        "foreign_keys": [],
    },
    {
        "name": "pathways",
        "path": "processed/nodes/pathways.tsv",
        "kind": "node",
        "node_type": "pathway",
        "expected_rows": 1_633,
        "fields": [
            field("canonical_id"),
            field("name"),
            field(
                "is_disease_pathway",
                "boolean",
            ),
        ],
        "primary_key": "canonical_id",
        "foreign_keys": [],
    },
    {
        "name": "tissues",
        "path": "processed/nodes/tissues.tsv",
        "kind": "node",
        "node_type": "tissue",
        "expected_rows": 301,
        "fields": [
            field("canonical_id"),
            field("name"),
        ],
        "primary_key": "canonical_id",
        "foreign_keys": [],
    },
    {
        "name": "gene_gene",
        "path": "processed/edges/gene_gene.tsv",
        "kind": "edge",
        "edge_type": "gene_gene",
        "expected_rows": 422_581,
        "fields": [
            field("source_id"),
            field("target_id"),
            field("relation"),
            field(
                "combined_score",
                "integer",
            ),
            field(
                "experimental_score",
                "integer",
            ),
            field(
                "database_score",
                "integer",
            ),
            field(
                "textmining_score",
                "integer",
            ),
            field(
                "supporting_protein_pairs",
                "integer",
            ),
        ],
        "primary_key": [
            "source_id",
            "target_id",
            "relation",
        ],
        "foreign_keys": [
            {
                "fields": "source_id",
                "reference": {
                    "resource": "genes",
                    "fields": "canonical_id",
                },
            },
            {
                "fields": "target_id",
                "reference": {
                    "resource": "genes",
                    "fields": "canonical_id",
                },
            },
        ],
    },
    {
        "name": "chemical_gene",
        "path": "processed/edges/chemical_gene.tsv",
        "kind": "edge",
        "edge_type": "chemical_gene",
        "expected_rows": 595_993,
        "fields": [
            field("source_id"),
            field("target_id"),
            field("relation"),
            field(
                "supporting_curated_rows",
                "integer",
            ),
            field(
                "supporting_pubmed_count",
                "integer",
            ),
            field("pubmed_ids"),
            field("interaction_actions"),
            field("gene_forms"),
        ],
        "primary_key": [
            "source_id",
            "target_id",
            "relation",
        ],
        "foreign_keys": [
            {
                "fields": "source_id",
                "reference": {
                    "resource": "chemicals",
                    "fields": "canonical_id",
                },
            },
            {
                "fields": "target_id",
                "reference": {
                    "resource": "genes",
                    "fields": "canonical_id",
                },
            },
        ],
    },
    {
        "name": "gene_pathway",
        "path": "processed/edges/gene_pathway.tsv",
        "kind": "edge",
        "edge_type": "gene_pathway",
        "expected_rows": 47_326,
        "fields": [
            field("source_id"),
            field("target_id"),
            field("relation"),
            field("evidence_codes"),
            field(
                "supporting_mapping_rows",
                "integer",
            ),
        ],
        "primary_key": [
            "source_id",
            "target_id",
            "relation",
        ],
        "foreign_keys": [
            {
                "fields": "source_id",
                "reference": {
                    "resource": "genes",
                    "fields": "canonical_id",
                },
            },
            {
                "fields": "target_id",
                "reference": {
                    "resource": "pathways",
                    "fields": "canonical_id",
                },
            },
        ],
    },
    {
        "name": "gene_tissue",
        "path": "processed/edges/gene_tissue.tsv",
        "kind": "edge",
        "edge_type": "gene_tissue",
        "expected_rows": 548_568,
        "fields": [
            field("source_id"),
            field("target_id"),
            field("relation"),
            field("bgee_gene_id"),
            field("call_quality"),
            field(
                "fdr",
                "number",
            ),
            field(
                "expression_score",
                "number",
            ),
            field(
                "expression_rank",
                "number",
            ),
            field(
                "self_observation_count",
                "integer",
            ),
            field(
                "descendant_observation_count",
                "integer",
            ),
        ],
        "primary_key": [
            "source_id",
            "target_id",
            "relation",
        ],
        "foreign_keys": [
            {
                "fields": "source_id",
                "reference": {
                    "resource": "genes",
                    "fields": "canonical_id",
                },
            },
            {
                "fields": "target_id",
                "reference": {
                    "resource": "tissues",
                    "fields": "canonical_id",
                },
            },
        ],
    },
]


def build_resource(
    spec: dict,
) -> dict:
    path = (
        DATA_ROOT
        / spec["path"]
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing release resource: {path}"
        )

    (
        observed_header,
        observed_rows,
    ) = inspect_tsv(
        path
    )

    expected_header = [
        item["name"]
        for item
        in spec["fields"]
    ]

    if (
        observed_header
        != expected_header
    ):
        raise ValueError(
            "Header mismatch for "
            f"{spec['name']}:\n"
            f"expected={expected_header}\n"
            f"observed={observed_header}"
        )

    if (
        observed_rows
        != spec[
            "expected_rows"
        ]
    ):
        raise ValueError(
            "Row-count mismatch for "
            f"{spec['name']}: "
            f"expected "
            f"{spec['expected_rows']:,}, "
            f"observed "
            f"{observed_rows:,}"
        )

    schema = {
        "fields":
            spec[
                "fields"
            ],
        "missingValues": [
            "NA"
        ],
        "primaryKey":
            spec[
                "primary_key"
            ],
    }

    if spec[
        "foreign_keys"
    ]:
        schema[
            "foreignKeys"
        ] = spec[
            "foreign_keys"
        ]

    resource = {
        "profile":
            "tabular-data-resource",
        "name":
            spec[
                "name"
            ],
        "path":
            spec[
                "path"
            ],
        "format":
            "tsv",
        "mediatype":
            "text/tab-separated-values",
        "encoding":
            "utf-8",
        "bytes":
            path.stat().st_size,
        "hash":
            (
                "sha256:"
                + sha256_file(
                    path
                )
            ),
        "rowCount":
            observed_rows,
        "dialect": {
            "delimiter":
                "\t",
            "lineTerminator":
                "\n",
            "quoteChar":
                "\"",
            "doubleQuote":
                True,
            "header":
                True,
        },
        "schema":
            schema,
        "netem": {
            "resource_kind":
                spec[
                    "kind"
                ],
        },
    }

    if (
        spec[
            "kind"
        ]
        == "node"
    ):
        resource[
            "netem"
        ][
            "node_type"
        ] = spec[
            "node_type"
        ]

    else:
        resource[
            "netem"
        ][
            "edge_type"
        ] = spec[
            "edge_type"
        ]

    return resource


def main() -> None:
    dataset = load_toml(
        DATASET_CONFIG
    )

    contract = load_toml(
        CONTRACT_CONFIG
    )

    sources_config = load_toml(
        SOURCES_CONFIG
    )

    core = dataset[
        "core_graph_v1"
    ]

    if (
        core[
            "status"
        ]
        != "frozen"
    ):
        raise ValueError(
            "Core Graph V1 must be "
            "frozen before release"
        )

    if (
        str(
            core[
                "version"
            ]
        )
        != "1.0.0"
    ):
        raise ValueError(
            "Unexpected release version"
        )

    resources = [
        build_resource(
            spec
        )
        for spec
        in RESOURCE_SPECS
    ]

    node_rows = sum(
        resource[
            "rowCount"
        ]
        for resource
        in resources
        if (
            resource[
                "netem"
            ][
                "resource_kind"
            ]
            == "node"
        )
    )

    edge_rows = sum(
        resource[
            "rowCount"
        ]
        for resource
        in resources
        if (
            resource[
                "netem"
            ][
                "resource_kind"
            ]
            == "edge"
        )
    )

    if (
        node_rows
        != EXPECTED_TOTAL_NODES
    ):
        raise ValueError(
            "Release node total mismatch: "
            f"{node_rows:,}"
        )

    if (
        edge_rows
        != EXPECTED_TOTAL_EDGES
    ):
        raise ValueError(
            "Release edge total mismatch: "
            f"{edge_rows:,}"
        )

    if (
        node_rows
        != int(
            core[
                "node_counts"
            ][
                "total"
            ]
        )
    ):
        raise ValueError(
            "Manifest node total does "
            "not match dataset config"
        )

    if (
        edge_rows
        != int(
            core[
                "edge_counts"
            ][
                "total"
            ]
        )
    ):
        raise ValueError(
            "Manifest edge total does "
            "not match dataset config"
        )

    config_files = {}

    for path in (
        DATASET_CONFIG,
        CONTRACT_CONFIG,
        SOURCES_CONFIG,
    ):
        config_files[
            path.name
        ] = {
            "sha256":
                sha256_file(
                    path
                ),
        }

    manifest = {
        "profile":
            "tabular-data-package",

        "name":
            "netem-rebuild",

        "title":
            "NetEM-Rebuild Core Graph V1",

        "description":
            (
                "Reproducible heterogeneous "
                "human biological network "
                "containing Gene, Chemical, "
                "Pathway, and Tissue nodes."
            ),

        "version":
            str(
                core[
                    "version"
                ]
            ),

        "keywords": [
            "heterogeneous-network",
            "disease-module",
            "gene",
            "chemical",
            "pathway",
            "tissue",
            "graph-neural-network",
        ],

        "sources": [
            {
                "title": "HGNC",
                "path":
                    "https://www.genenames.org/",
            },
            {
                "title": "STRING",
                "path":
                    "https://string-db.org/",
            },
            {
                "title":
                    "Comparative Toxicogenomics Database",
                "path":
                    "https://ctdbase.org/",
            },
            {
                "title": "Reactome",
                "path":
                    "https://reactome.org/",
            },
            {
                "title": "Bgee",
                "path":
                    "https://www.bgee.org/",
            },
        ],

        "resources":
            resources,

        "netem": {
            "organism": {
                "name":
                    dataset[
                        "dataset"
                    ][
                        "organism"
                    ],
                "taxon_id":
                    dataset[
                        "dataset"
                    ][
                        "taxon_id"
                    ],
            },

            "core_graph": {
                "status":
                    core[
                        "status"
                    ],
                "version":
                    str(
                        core[
                            "version"
                        ]
                    ),
                "node_counts":
                    json_safe(
                        core[
                            "node_counts"
                        ]
                    ),
                "edge_counts":
                    json_safe(
                        core[
                            "edge_counts"
                        ]
                    ),
                "coverage":
                    json_safe(
                        core[
                            "coverage"
                        ]
                    ),
                "connectivity":
                    json_safe(
                        core[
                            "connectivity"
                        ]
                    ),
                "annotations":
                    json_safe(
                        core[
                            "annotations"
                        ]
                    ),
            },

            "data_contract": {
                "name":
                    contract[
                        "contract"
                    ][
                        "name"
                    ],
                "version":
                    str(
                        contract[
                            "contract"
                        ][
                            "version"
                        ]
                    ),
            },

            "config_files":
                config_files,

            "source_snapshots":
                json_safe(
                    sources_config[
                        "sources"
                    ]
                ),
        },
    }

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    serialized = (
        json.dumps(
            manifest,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    OUTPUT.write_text(
        serialized,
        encoding="utf-8",
    )

    # ----------------------------------------------
    # Round-trip validation
    # ----------------------------------------------

    loaded = json.loads(
        OUTPUT.read_text(
            encoding="utf-8"
        )
    )

    if (
        loaded
        != manifest
    ):
        raise ValueError(
            "Manifest JSON round-trip "
            "validation failed"
        )

    if (
        len(
            loaded[
                "resources"
            ]
        )
        != 8
    ):
        raise ValueError(
            "Manifest must contain "
            "exactly 8 resources"
        )

    print(
        "Release resources verified:",
        len(resources),
    )

    print(
        "Node rows:",
        f"{node_rows:,}",
    )

    print(
        "Edge rows:",
        f"{edge_rows:,}",
    )

    print(
        "Dataset version:",
        core[
            "version"
        ],
    )

    print(
        "Data contract version:",
        contract[
            "contract"
        ][
            "version"
        ],
    )

    print(
        "Wrote:",
        OUTPUT.relative_to(
            PROJECT_ROOT
        ),
    )

    print(
        "Manifest SHA-256:",
        sha256_file(
            OUTPUT
        ),
    )

    print()

    print(
        "Core Graph V1 release manifest: OK"
    )


if __name__ == "__main__":
    main()