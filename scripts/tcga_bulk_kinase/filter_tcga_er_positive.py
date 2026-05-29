#!/usr/bin/env python3
"""Filter a TCGA BRCA expression matrix to ER-positive primary-tumor samples.

Inputs
------
1. A tab-delimited expression matrix with genes as rows and TCGA sample barcodes
   as columns.
2. A CSV containing ER-positive TCGA case IDs in a CLID column.

The ER-positive file is expected to contain case-level barcodes without sample
suffixes, for example TCGA-XX-YYYY. This script appends a configurable suffix
(default: -01) to retain primary tumor samples from the expression matrix.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/path/to/TCGA_data"),
        help="Directory containing the TCGA input files and outputs.",
    )
    parser.add_argument(
        "--er-positive-csv",
        type=Path,
        default=None,
        help="CSV with a CLID column listing ER-positive TCGA case IDs.",
    )
    parser.add_argument(
        "--expression-file",
        type=Path,
        default=None,
        help="Tab-delimited expression matrix; genes are rows and samples are columns.",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=None,
        help="Filtered expression matrix to write.",
    )
    parser.add_argument(
        "--clid-column",
        default="CLID",
        help="Column in the ER-positive CSV containing TCGA case IDs.",
    )
    parser.add_argument(
        "--sample-suffix",
        default="-01",
        help="Sample suffix to append to case IDs. Default keeps primary tumors.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir
    er_positive_csv = args.er_positive_csv or base_dir / "ER_possitive_IHC_TCGA_patients.csv"
    expression_file = args.expression_file or base_dir / "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
    output_file = args.output_file or base_dir / "filtered_gene_expression_data.txt"

    er_positive_patients = pd.read_csv(er_positive_csv)
    if args.clid_column not in er_positive_patients.columns:
        raise ValueError(
            f"Column {args.clid_column!r} was not found in {er_positive_csv}. "
            f"Available columns: {list(er_positive_patients.columns)}"
        )

    clids = (
        er_positive_patients[args.clid_column]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
    )
    clids_to_keep = [f"{clid}{args.sample_suffix}" for clid in clids]

    expression_df = pd.read_csv(expression_file, sep="\t", index_col=0)
    expression_df.index = expression_df.index.map(str)

    matching_clids = [clid for clid in clids_to_keep if clid in expression_df.columns]
    missing_clids = sorted(set(clids_to_keep) - set(matching_clids))
    if not matching_clids:
        raise ValueError(
            "No ER-positive sample IDs matched the expression matrix columns. "
            "Check the CLID format and sample suffix."
        )

    filtered_df = expression_df.loc[:, matching_clids]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    filtered_df.to_csv(output_file, sep="\t")

    print(f"Input ER-positive cases: {len(clids)}")
    print(f"Matched expression columns: {len(matching_clids)}")
    print(f"Missing requested columns: {len(missing_clids)}")
    print(f"Wrote: {output_file}")


if __name__ == "__main__":
    main()
