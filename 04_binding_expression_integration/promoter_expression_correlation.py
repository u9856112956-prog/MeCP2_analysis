#!/usr/bin/env python3
"""
Promoter-Specific MeCP2 Binding vs Absolute Expression Changes Analysis (DEA genes only)

This script quantifies how differences in MeCP2 binding between exogenous (Exo) and
endogenous (Endo) samples correlate with gene expression changes. It tests whether
promoter binding differences show stronger correlation with Exo vs Endo expression
log2 fold change than gene body regions.

This script performs the same analysis as promoter_expression_absolute_correlation.py
but restricts the analysis to only differentially expressed genes defined
in the DEA filtered files.

Key features:
- Uses absolute values of log2FoldChange to focus on dysregulation magnitude
- Restricts analysis to only DEA genes for more focused hypothesis testing
- Analyzes how MeCP2 binding changes predict expression change magnitude

Usage example:
    python promoter_expression_absolute_correlation_deaonly.py \
        --integrated-file /path/to/integrated_gene_body_rnaseq_NSC.csv \
        --dea-filtered-file /path/to/DEA_NSC_filtered.csv \
        --output-dir /path/to/output \
        --cell-type NSC
"""

import argparse
import logging
import math
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

EPS = 1e-6

REGION_LABELS: Dict[str, str] = {
    "promoter": "Promoter",
    "gene_body_5prime": "5' Gene Body",
    "gene_body_middle": "Middle Gene Body",
    "gene_body_3prime": "3' Gene Body",
    "gene_body_full": "Full Gene Body",
}

SUBSET_DEFINITIONS = {
    "all_dea": lambda df: df,  # All DEA genes
    "upregulated": lambda df: df[df["log2FoldChange"] > 0],
    "downregulated": lambda df: df[df["log2FoldChange"] < 0],
}


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def load_dea_filtered_genes(path: Path) -> pd.DataFrame:
    """Load the DEA filtered genes file."""
    logging.info("Loading DEA filtered genes from %s", path)
    df = pd.read_csv(path)

    # Ensure gene column is clean (remove version numbers)
    df["gene_clean"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    logging.info("Loaded %d DEA filtered genes", len(df))
    return df


def filter_integrated_by_dea(integrated_df: pd.DataFrame, dea_df: pd.DataFrame) -> pd.DataFrame:
    """Filter integrated dataset to only include DEA genes."""
    logging.info("Filtering integrated dataset to include only DEA genes")

    # Get list of DEA gene identifiers
    dea_genes = set(dea_df["gene_clean"].tolist())

    # Filter integrated data by DEA genes
    # Try different identifier columns
    if "gene" in integrated_df.columns:
        gene_col = "gene"
    elif "gene_identifier" in integrated_df.columns:
        gene_col = "gene_identifier"
    elif "gene_name" in integrated_df.columns:
        gene_col = "gene_name"
    else:
        raise ValueError("Could not find suitable gene identifier column in integrated data")

    # Clean gene identifiers in integrated data
    integrated_df = integrated_df.copy()
    integrated_df["gene_clean"] = integrated_df[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)

    # Filter to only DEA genes
    filtered_df = integrated_df[integrated_df["gene_clean"].isin(dea_genes)].copy()

    logging.info("Filtered from %d to %d genes (DEA only)", len(integrated_df), len(filtered_df))

    if len(filtered_df) == 0:
        raise ValueError("No overlapping genes found between integrated data and DEA filtered genes")

    return filtered_df


def _compute_corr(
    x: pd.Series,
    y: pd.Series,
    method: str,
) -> Tuple[float, float]:
    valid = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    x_valid = x[valid]
    y_valid = y[valid]

    if len(x_valid) < 3 or x_valid.nunique() < 2 or y_valid.nunique() < 2:
        return np.nan, np.nan

    if method == "pearson":
        corr, pval = stats.pearsonr(x_valid, y_valid)
    else:
        corr, pval = stats.spearmanr(x_valid, y_valid)
    return corr, pval


def _compute_regression(x: pd.Series, y: pd.Series) -> Tuple[float, float, float]:
    valid = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    x_valid = x[valid]
    y_valid = y[valid]

    if len(x_valid) < 3 or x_valid.nunique() < 2 or y_valid.nunique() < 2:
        return np.nan, np.nan, np.nan

    lr = stats.linregress(x_valid, y_valid)
    return lr.slope, lr.pvalue, lr.rvalue ** 2


def _fisher_ci(corr: float, n: int, alpha: float) -> Tuple[float, float]:
    if not np.isfinite(corr) or n <= 3:
        return np.nan, np.nan

    corr = max(min(corr, 0.999999), -0.999999)
    z = np.arctanh(corr)
    se = 1 / math.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    lower = np.tanh(z - z_crit * se)
    upper = np.tanh(z + z_crit * se)
    return lower, upper


def compute_region_correlations(
    integrated_df: pd.DataFrame,
    method: str,
    alpha: float,
    min_genes: int,
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    records = []
    scatter_data: Dict[str, pd.DataFrame] = {}

    for region, label in REGION_LABELS.items():
        exo_col = f"{region}_exo_signal"
        endo_col = f"{region}_endo_signal"

        if exo_col not in integrated_df.columns or endo_col not in integrated_df.columns:
            logging.warning("Skipping region %s: required columns missing", region)
            continue

        region_df = integrated_df[[exo_col, endo_col, "log2FoldChange", "padj"]].copy()
        region_df["log2_ratio"] = np.log2(
            (region_df[exo_col] + EPS) / (region_df[endo_col] + EPS)
        )
        # Key change: use absolute value of log2FoldChange
        region_df["abs_log2FoldChange"] = region_df["log2FoldChange"].abs()
        region_df["delta_signal"] = region_df[endo_col] - region_df[exo_col]
        region_df = region_df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["log2FoldChange", "log2_ratio"]
        )

        if region_df.empty:
            logging.warning("Skipping region %s: no valid data after filtering", region)
            continue

        # Store both absolute and original values for scatter plots
        scatter_data[region] = region_df[["log2_ratio", "abs_log2FoldChange", "log2FoldChange"]].copy()

        for subset_name, subset_fn in SUBSET_DEFINITIONS.items():
            subset = subset_fn(region_df).copy()
            subset = subset.dropna(subset=["log2_ratio", "abs_log2FoldChange"])

            if len(subset) < min_genes:
                logging.debug(
                    "Skipping region %s subset %s: only %d genes (min=%d)",
                    region,
                    subset_name,
                    len(subset),
                    min_genes,
                )
                continue

            # Correlate with absolute log2FoldChange
            corr_ratio, p_ratio = _compute_corr(
                subset["log2_ratio"], subset["abs_log2FoldChange"], method
            )
            corr_delta, p_delta = _compute_corr(
                subset["delta_signal"], subset["abs_log2FoldChange"], method
            )
            slope, slope_p, r_sq = _compute_regression(
                subset["log2_ratio"], subset["abs_log2FoldChange"]
            )
            ci_lower, ci_upper = _fisher_ci(corr_ratio, len(subset), alpha)

            # Also compute correlation with signed log2FC for comparison
            corr_signed, p_signed = _compute_corr(
                subset["log2_ratio"], subset["log2FoldChange"], method
            )

            records.append(
                {
                    "region": region,
                    "region_label": label,
                    "subset": subset_name,
                    "n_genes": len(subset),
                    "correlation_method": method,
                    "correlation_log2_ratio_abs": corr_ratio,
                    "pvalue_log2_ratio_abs": p_ratio,
                    "correlation_log2_ratio_signed": corr_signed,
                    "pvalue_log2_ratio_signed": p_signed,
                    "correlation_delta_abs": corr_delta,
                    "pvalue_delta_abs": p_delta,
                    "fisher_ci_lower": ci_lower,
                    "fisher_ci_upper": ci_upper,
                    "linear_slope": slope,
                    "linear_pvalue": slope_p,
                    "r_squared": r_sq,
                    "mean_log2_ratio": subset["log2_ratio"].mean(),
                    "mean_delta_signal": subset["delta_signal"].mean(),
                    "mean_abs_expression_log2FC": subset["abs_log2FoldChange"].mean(),
                    "mean_signed_expression_log2FC": subset["log2FoldChange"].mean(),
                }
            )

    results_df = pd.DataFrame(records)
    return results_df, scatter_data


def create_correlation_plot(
    results_df: pd.DataFrame,
    alpha: float,
    output_dir: Path,
    suffix: str,
) -> None:
    subset_df = results_df[results_df["subset"] == "all_dea"].copy()
    if subset_df.empty:
        logging.warning("No 'all_dea' subset correlations available; skipping summary plot.")
        return

    subset_df = subset_df.sort_values("correlation_log2_ratio_abs", ascending=False)
    x_positions = np.arange(len(subset_df))
    lower_err = (
        subset_df["correlation_log2_ratio_abs"] - subset_df["fisher_ci_lower"]
    ).clip(lower=0).fillna(0.0)
    upper_err = (
        subset_df["fisher_ci_upper"] - subset_df["correlation_log2_ratio_abs"]
    ).clip(lower=0).fillna(0.0)
    yerr = np.vstack([lower_err.to_numpy(), upper_err.to_numpy()])

    colors = [
        "#d62728" if region == "promoter" else "#1f77b4"
        for region in subset_df["region"]
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(
        x_positions,
        subset_df["correlation_log2_ratio_abs"],
        yerr=yerr,
        fmt="o",
        ecolor="gray",
        elinewidth=2,
        capsize=4,
        mec="black",
        mfc="white",
    )
    ax.scatter(
        x_positions,
        subset_df["correlation_log2_ratio_abs"],
        s=120,
        color=colors,
        edgecolor="black",
        zorder=3,
    )
    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.6)

    for idx, row in enumerate(subset_df.itertuples()):
        if np.isfinite(row.pvalue_log2_ratio_abs) and row.pvalue_log2_ratio_abs < alpha:
            ax.text(
                idx,
                row.correlation_log2_ratio_abs + 0.03,
                "*",
                ha="center",
                va="bottom",
                color=colors[idx],
                fontsize=14,
                fontweight="bold",
            )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(subset_df["region_label"], rotation=45, ha="right")
    ax.set_ylabel(
        f"{subset_df['correlation_method'].iloc[0].title()} correlation "
        "(log2 MeCP2 Exo/Endo vs |log2 Expression Exo/Endo|)"
    )
    ax.set_title("Promoter vs Gene Body Correlation with Absolute Expression Changes (DEA genes only)")
    ax.set_ylim(
        min(-0.1, subset_df["correlation_log2_ratio_abs"].min() - 0.1),
        max(0.5, subset_df["correlation_log2_ratio_abs"].max() + 0.1),
    )
    plt.tight_layout()

    fig_path = output_dir / f"binding_absolute_expression_correlation_overview_deaonly_{suffix}.png"
    plt.savefig(fig_path, dpi=300)
    plt.savefig(fig_path.with_suffix(".pdf"))
    plt.close(fig)


def create_scatter_grid(
    scatter_data: Dict[str, pd.DataFrame],
    output_dir: Path,
    suffix: str,
    sample_size: int = 5000,
    log2fc_outlier_threshold: float = 15.0,
) -> None:
    if not scatter_data:
        logging.warning("No scatter data available; skipping scatter grid.")
        return

    regions_in_plot = [region for region in REGION_LABELS if region in scatter_data]
    n_regions = len(regions_in_plot)
    ncols = 3
    nrows = math.ceil(n_regions / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes_flat = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

    for idx, region in enumerate(regions_in_plot):
        ax = axes_flat[idx]
        data = scatter_data[region]
        if data.empty:
            ax.set_visible(False)
            continue

        if len(data) > sample_size:
            data_plot = data.sample(sample_size, random_state=42)
        else:
            data_plot = data

        # Filter outliers: remove extreme values from both axes
        # Filter y-axis (abs_log2FoldChange)
        y_outlier_filter = data_plot['abs_log2FoldChange'] <= log2fc_outlier_threshold
        # Filter x-axis (log2_ratio of exo/endo signals) - use same threshold for consistency
        x_outlier_filter = (data_plot['log2_ratio'] <= log2fc_outlier_threshold) & (data_plot['log2_ratio'] >= -log2fc_outlier_threshold)
        
        # Apply both filters
        combined_filter = y_outlier_filter & x_outlier_filter
        data_filtered = data_plot[combined_filter]
        
        n_outliers = len(data_plot) - len(data_filtered)
        if n_outliers > 0:
            logging.info(f"Filtered {n_outliers} outliers with |log2FC| > {log2fc_outlier_threshold} or |log2 ratio| > {log2fc_outlier_threshold} for {region} region")

        sns.regplot(
            data=data_filtered,
            x="log2_ratio",
            y="abs_log2FoldChange",
            scatter_kws={"s": 16, "alpha": 0.35, "edgecolor": "none"},
            line_kws={"color": "#d62728"},
            ax=ax,
        )
        ax.set_xlabel("log2(MeCP2 Exo / Endo) signal")
        ax.set_ylabel("|log2(Expression Exo / Endo)|")
        ax.set_title(f"{REGION_LABELS[region]} (DEA genes, outliers removed)")

        ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)

    for idx in range(n_regions, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()
    fig_path = output_dir / f"binding_absolute_expression_scatter_grid_deaonly_{suffix}.png"
    plt.savefig(fig_path, dpi=300)
    plt.savefig(fig_path.with_suffix(".pdf"))
    plt.close(fig)


def write_summary(
    results_df: pd.DataFrame,
    alpha: float,
    output_dir: Path,
    suffix: str,
    n_dea_genes: int,
) -> None:
    summary_path = output_dir / f"binding_absolute_expression_summary_deaonly_{suffix}.txt"
    subset_df = results_df[results_df["subset"] == "all_dea"].copy()

    with open(summary_path, "w") as handle:
        handle.write("Promoter vs Absolute Gene Expression Changes Correlation Summary (DEA genes only)\n")
        handle.write("=" * 80 + "\n\n")
        handle.write("Analysis focus: Correlation with magnitude of gene dysregulation in DEA genes\n")
        handle.write("(using absolute values of log2 fold changes)\n\n")
        handle.write(f"Analysis restricted to {n_dea_genes} differentially expressed genes\n")
        handle.write(f"Significance threshold (alpha): {alpha}\n")
        method = (
            subset_df['correlation_method'].iloc[0]
            if not subset_df.empty
            else 'NA'
        )
        handle.write(f"Correlation method: {method}\n\n")

        if subset_df.empty:
            handle.write("No correlation results available for DEA genes.\n")
            return

        promoter_row = subset_df[subset_df["region"] == "promoter"]
        if not promoter_row.empty:
            row = promoter_row.iloc[0]
            handle.write("Promoter correlation with absolute expression changes (DEA genes):\n")
            handle.write(
                f"  r = {row.correlation_log2_ratio_abs:.3f}, p = {row.pvalue_log2_ratio_abs:.2e}, "
                f"n = {row.n_genes}, 95% CI = [{row.fisher_ci_lower:.3f}, {row.fisher_ci_upper:.3f}]\n"
            )
            handle.write(
                f"  For comparison, signed correlation: r = {row.correlation_log2_ratio_signed:.3f}, "
                f"p = {row.pvalue_log2_ratio_signed:.2e}\n\n"
            )
        else:
            handle.write("Promoter correlation not available.\n\n")

        handle.write("Gene body regions (absolute expression changes):\n")
        grouped = subset_df[subset_df["region"] != "promoter"]
        if grouped.empty:
            handle.write("  No gene body correlation results available.\n")
        else:
            for row in grouped.itertuples():
                handle.write(
                    f"  {row.region_label}: r = {row.correlation_log2_ratio_abs:.3f}, "
                    f"p = {row.pvalue_log2_ratio_abs:.2e}, n = {row.n_genes}\n"
                )

        significant_non_promoter = grouped[
            grouped["pvalue_log2_ratio_abs"] < alpha
        ]["region_label"].tolist()
        if significant_non_promoter:
            handle.write(
                "\nGene body regions with significant correlations (absolute changes): "
                + ", ".join(significant_non_promoter)
                + "\n"
            )
        else:
            handle.write(
                "\nNo gene body regions show significant correlation with absolute expression changes.\n"
            )

        handle.write("\nExpression direction breakdown (DEA genes):\n")
        upregulated_df = results_df[results_df["subset"] == "upregulated"]
        downregulated_df = results_df[results_df["subset"] == "downregulated"]

        if not upregulated_df.empty:
            promoter_up = upregulated_df[upregulated_df["region"] == "promoter"]
            if not promoter_up.empty:
                row = promoter_up.iloc[0]
                handle.write(
                    f"  Upregulated genes (n={row.n_genes}): promoter |r| = {row.correlation_log2_ratio_abs:.3f}, "
                    f"p = {row.pvalue_log2_ratio_abs:.2e}\n"
                )

        if not downregulated_df.empty:
            promoter_down = downregulated_df[downregulated_df["region"] == "promoter"]
            if not promoter_down.empty:
                row = promoter_down.iloc[0]
                handle.write(
                    f"  Downregulated genes (n={row.n_genes}): promoter |r| = {row.correlation_log2_ratio_abs:.3f}, "
                    f"p = {row.pvalue_log2_ratio_abs:.2e}\n"
                )

        handle.write("\nInterpretation (DEA genes, absolute changes):\n")
        if not promoter_row.empty:
            if promoter_row["pvalue_log2_ratio_abs"].iloc[0] < alpha:
                handle.write(
                    "- Promoter binding differences are significantly associated with the MAGNITUDE "
                    "of expression changes in differentially expressed genes (regardless of direction).\n"
                )
            else:
                handle.write(
                    "- Promoter binding differences do not significantly correlate with expression magnitude "
                    "in DEA genes; may indicate other regulatory mechanisms control change magnitude.\n"
                )
        else:
            handle.write(
                "- Promoter results unavailable; ensure promoter columns were present in the input file.\n"
            )

        if significant_non_promoter:
            handle.write(
                "- At least one non-promoter region also correlates with expression magnitude; "
                "inspect detailed table for effect sizes.\n"
            )
        else:
            handle.write(
                "- Gene body regions lack significant correlations with expression magnitude in DEA genes, "
                "supporting promoter-specific effects on gene dysregulation magnitude.\n"
            )

        handle.write(
            "\nNote: This analysis focuses on dysregulation magnitude in already-changing genes. "
            "Higher correlations indicate that MeCP2 binding changes predict how much "
            "a differentially expressed gene's expression will change, regardless of direction.\n"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assess MeCP2 binding correlation with absolute expression changes in DEA genes only."
    )
    parser.add_argument(
        "--integrated-file",
        type=Path,
        required=True,
        help="CSV file produced by integrate_gene_body_with_rnaseq.py.",
    )
    parser.add_argument(
        "--dea-filtered-file",
        type=Path,
        required=True,
        help="CSV file containing differentially expressed genes (DEA_NSC_filtered.csv or DEA_NEU_filtered.csv).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to save correlation results and plots.",
    )
    parser.add_argument(
        "--cell-type",
        default=None,
        help="Optional cell type label for output filenames.",
    )
    parser.add_argument(
        "--correlation-method",
        choices=["pearson", "spearman"],
        default="spearman",
        help="Correlation metric to use.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold for interpreting correlations.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=50,
        help="Minimum number of genes required per subset to compute statistics.",
    )
    parser.add_argument(
        "--log2fc-outlier-threshold",
        type=float,
        default=15.0,
        help="Filter out genes with |log2FC| > threshold (default 15.0). Use to remove extreme outliers.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_logging(args.verbose)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load integrated data
    integrated_df = pd.read_csv(args.integrated_file)
    logging.info(
        "Loaded integrated dataset with %d genes from %s",
        len(integrated_df),
        args.integrated_file,
    )

    # Load DEA filtered genes
    dea_df = load_dea_filtered_genes(args.dea_filtered_file)

    # Filter integrated data to only include DEA genes
    filtered_df = filter_integrated_by_dea(integrated_df, dea_df)

    if "log2FoldChange" not in filtered_df.columns or "padj" not in filtered_df.columns:
        raise ValueError(
            "Filtered dataset must contain 'log2FoldChange' and 'padj' columns."
        )

    suffix = (args.cell_type or "overall").lower()

    results_df, scatter_data = compute_region_correlations(
        integrated_df=filtered_df,
        method=args.correlation_method,
        alpha=args.alpha,
        min_genes=args.min_genes,
    )

    if results_df.empty:
        logging.error("No correlation results were produced. Check input data and thresholds.")
        return

    results_path = output_dir / f"binding_absolute_expression_correlation_deaonly_{suffix}.csv"
    results_df.to_csv(results_path, index=False)
    logging.info("Saved correlation table to %s", results_path)

    legacy_path = output_dir / f"binding_expression_correlation_{suffix}.csv"
    results_df.to_csv(legacy_path, index=False)
    logging.info("Saved legacy-formatted correlation table to %s", legacy_path)

    create_correlation_plot(results_df, args.alpha, output_dir, suffix)
    create_scatter_grid(scatter_data, output_dir, suffix, log2fc_outlier_threshold=args.log2fc_outlier_threshold)
    write_summary(results_df, args.alpha, output_dir, suffix, len(dea_df))

    logging.info("Promoter absolute expression correlation analysis (DEA only) complete.")


if __name__ == "__main__":
    main()
