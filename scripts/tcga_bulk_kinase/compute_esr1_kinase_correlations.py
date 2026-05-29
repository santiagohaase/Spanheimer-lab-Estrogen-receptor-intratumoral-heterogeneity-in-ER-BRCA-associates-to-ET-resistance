#!/usr/bin/env python3
"""Correlate ESR1 expression with TCGA ER-positive kinase ssGSEA NES values."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/path/to/TCGA_data"),
        help="Directory containing TCGA input files and outputs.",
    )
    parser.add_argument("--expression-file", type=Path, default=None)
    parser.add_argument("--nes-matrix", type=Path, default=None)
    parser.add_argument("--output-file", type=Path, default=None)
    parser.add_argument("--gene", default="ESR1", help="Gene symbol to correlate with kinase NES values.")
    return parser.parse_args()


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Return Benjamini-Hochberg adjusted q-values, preserving the input index."""
    p = p_values.astype(float)
    valid = p.notna()
    q = pd.Series(np.nan, index=p.index, dtype=float)
    if valid.sum() == 0:
        return q

    ranked = p[valid].sort_values()
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    q.loc[adjusted.index] = adjusted
    return q


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir
    expression_file = args.expression_file or base_dir / "filtered_gene_expression_data.txt"
    nes_matrix_file = args.nes_matrix or base_dir / "ssGSEA_results" / "ssgsea_NES_matrix.tsv"
    output_file = args.output_file or base_dir / "ssGSEA_results" / "ESR1_kinase_signature_correlations.tsv"

    expression_df = pd.read_csv(expression_file, sep="\t", index_col=0)
    expression_df.index = expression_df.index.map(str)
    if args.gene not in expression_df.index:
        raise ValueError(f"{args.gene!r} was not found in the expression matrix index.")

    gene_expr = pd.to_numeric(expression_df.loc[args.gene], errors="coerce")
    signature_scores = pd.read_csv(nes_matrix_file, sep="\t", index_col=0)
    signature_scores = signature_scores.apply(pd.to_numeric, errors="coerce")

    common_samples = signature_scores.index.intersection(gene_expr.index)
    if len(common_samples) < 3:
        raise ValueError(
            f"Only {len(common_samples)} shared samples were found between expression and NES matrices."
        )

    gene_expr = gene_expr.loc[common_samples]
    signature_scores = signature_scores.loc[common_samples]

    rows = []
    for signature in signature_scores.columns:
        pair = pd.concat(
            [gene_expr.rename(args.gene), signature_scores[signature].rename("NES")],
            axis=1,
        ).dropna()
        if pair.shape[0] < 3 or pair[args.gene].nunique() < 2 or pair["NES"].nunique() < 2:
            r, p = np.nan, np.nan
        else:
            r, p = pearsonr(pair[args.gene], pair["NES"])
        rows.append(
            {
                "Signature": signature,
                "Pearson_r": r,
                "p_value": p,
                "n_samples": pair.shape[0],
            }
        )

    results_df = pd.DataFrame(rows)
    results_df["BH_q_value"] = benjamini_hochberg(results_df["p_value"])
    results_df = results_df.sort_values("Pearson_r", ascending=True, na_position="last")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_file, sep="\t", index=False)

    print(f"Shared samples used: {len(common_samples)}")
    print(f"Wrote correlations: {output_file}")
    print(results_df.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
