from scipy.stats import friedmanchisquare
import pandas as pd


def load_scores(csv_file):
    """
    Load the alignment scores CSV.

    Parameters
    ----------
    csv_file : str
        Path to alignment_scores.csv

    Returns
    -------
    pandas.DataFrame
    """
    return pd.read_csv(csv_file)


def prepare_metric_data(df, metric):
    """
    Convert the long-format dataframe into one row per dataset.

    Parameters
    ----------
    df : pandas.DataFrame
    metric : str
        "SP" or "TC"

    Returns
    -------
    pandas.DataFrame
        Index = Dataset
        Columns = aligners
    """

    if metric not in ["SP", "TC"]:
        raise ValueError("Metric must be 'SP' or 'TC'.")

    pivot = df.pivot(
        index="Dataset",
        columns="Aligner",
        values=metric
    )

    return pivot.sort_index()


def run_friedman_test(metric_table):
    """
    Perform the Friedman test.

    Parameters
    ----------
    metric_table : pandas.DataFrame
        Output from prepare_metric_data()

    Returns
    -------
    statistic : float
    pvalue : float
    """

    statistic, pvalue = friedmanchisquare(
        metric_table["clustal"],
        metric_table["mafft"],
        metric_table["muscle"]
    )

    return statistic, pvalue


def kendalls_w(statistic, n_datasets, n_aligners):
    """
    Compute Kendall's coefficient of concordance (W).

    W = chi^2 / (N * (k - 1))
    """

    return statistic / (n_datasets * (n_aligners - 1))


def analyze_metric(csv_file, metric):
    """
    Complete Friedman analysis for a metric.

    Parameters
    ----------
    csv_file : str
        Path to alignment_scores.csv
    metric : str
        "SP" or "TC"

    Returns
    -------
    dict containing:
        metric
        friedman_chi2
        pvalue
        kendalls_w
    """

    df = load_scores(csv_file)

    metric_table = prepare_metric_data(df, metric)

    statistic, pvalue = run_friedman_test(metric_table)

    w = kendalls_w(
        statistic,
        n_datasets=len(metric_table),
        n_aligners=metric_table.shape[1]
    )

    return {
        "metric": metric,
        "friedman_chi2": statistic,
        "pvalue": pvalue,
        "kendalls_w": w
    }