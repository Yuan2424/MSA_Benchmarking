from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


def perform_posthoc(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """
    Perform pairwise Wilcoxon signed-rank tests for a given metric.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns:
            Dataset, Aligner, SP, TC

    metric : str
        Either "SP" or "TC".

    Returns
    -------
    pandas.DataFrame
        Pairwise Wilcoxon test results with Holm-adjusted p-values.
    """

    metric = metric.upper()

    if metric not in ("SP", "TC"):
        raise ValueError("Metric must be either 'SP' or 'TC'.")

    # Convert long-format dataframe into paired columns:
    #
    # Dataset     clustal     mafft     muscle
    # BB11035      ...         ...        ...
    #
    pivot = df.pivot(
        index="Dataset",
        columns="Aligner",
        values=metric
    ).sort_index()

    aligners = sorted(pivot.columns)

    results = []
    raw_p_values = []

    # Compare every pair of aligners
    for aligner1, aligner2 in combinations(aligners, 2):

        statistic, p_value = wilcoxon(
            pivot[aligner1],
            pivot[aligner2],
            alternative="two-sided"
        )

        results.append({
            "Metric": metric,
            "Comparison": f"{aligner1} vs {aligner2}",
            "Statistic": statistic,
            "Raw p-value": p_value
        })

        raw_p_values.append(p_value)

    # Holm correction for multiple comparisons
    reject, corrected_p, _, _ = multipletests(
        raw_p_values,
        alpha=0.05,
        method="holm"
    )

    for i in range(len(results)):
        results[i]["Adjusted p-value"] = corrected_p[i]
        results[i]["Significant"] = reject[i]

    return pd.DataFrame(results)


def run_posthoc(
    input_csv: str = "Results/alignment_scores.csv",
    output_csv: str = "Results/wilcoxon.csv"
):
    """
    Run post-hoc Wilcoxon tests for both SP and TC metrics and save
    the results.

    Parameters
    ----------
    input_csv : str
        Path to alignment_scores.csv

    output_csv : str
        Path where Wilcoxon results should be written.
    """

    df = pd.read_csv(input_csv)

    sp_results = perform_posthoc(df, "SP")
    tc_results = perform_posthoc(df, "TC")

    results = pd.concat(
        [sp_results, tc_results],
        ignore_index=True
    )

    # Ensure Results directory exists
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)

    results.to_csv(output_csv, index=False)

    print(f"Wilcoxon results saved to: {output_csv}")