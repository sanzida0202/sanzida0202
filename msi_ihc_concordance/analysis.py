"""Core analysis helpers for the MSI vs. IHC concordance project."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm


@dataclass
class ConcordanceSummary:
    """Container for key concordance metrics."""

    total_cases: int
    concordant_cases: int
    discordant_cases: int
    concordance_rate: float
    discordance_rate: float
    contingency_table: Dict[str, Dict[str, int]]
    discordant_case_ids: Iterable[str]

    def to_dict(self) -> Dict[str, object]:
        """Convert the summary to a JSON-serialisable dictionary."""

        return {
            "total_cases": self.total_cases,
            "concordant_cases": self.concordant_cases,
            "discordant_cases": self.discordant_cases,
            "concordance_rate": self.concordance_rate,
            "discordance_rate": self.discordance_rate,
            "contingency_table": self.contingency_table,
            "discordant_case_ids": list(self.discordant_case_ids),
        }


def load_msi_calls(path: Path | str) -> pd.DataFrame:
    """Load MSI calls from a CSV file."""

    msi_df = pd.read_csv(path)
    msi_df.columns = [col.strip().lower() for col in msi_df.columns]
    expected_columns = {"sample_id", "msi_score", "msi_call"}
    missing = expected_columns.difference(msi_df.columns)
    if missing:
        raise ValueError(f"MSI call file is missing columns: {sorted(missing)}")

    msi_df["msi_call"] = msi_df["msi_call"].str.upper().str.strip()
    msi_df["msi_score"] = pd.to_numeric(msi_df["msi_score"], errors="coerce")
    if msi_df["msi_score"].isna().any():
        raise ValueError("All MSI scores must be numeric")

    return msi_df


def load_ihc_results(path: Path | str) -> pd.DataFrame:
    """Load IHC staining results from a CSV file."""

    ihc_df = pd.read_csv(path)
    ihc_df.columns = [col.strip() for col in ihc_df.columns]
    required_markers = {"MLH1", "MSH2", "MSH6", "PMS2"}
    missing_markers = required_markers.difference(set(ihc_df.columns))
    if missing_markers:
        raise ValueError(f"IHC results are missing markers: {sorted(missing_markers)}")

    return ihc_df


def _classify_ihc(row: pd.Series) -> Tuple[str, Tuple[str, ...]]:
    markers = ("MLH1", "MSH2", "MSH6", "PMS2")
    lost = tuple(marker for marker in markers if str(row[marker]).strip().lower() == "lost")
    status = "dMMR" if lost else "pMMR"
    return status, lost


def merge_assay_results(msi_df: pd.DataFrame, ihc_df: pd.DataFrame) -> pd.DataFrame:
    """Combine MSI calls with IHC staining results."""

    merged = pd.merge(msi_df, ihc_df, on="sample_id", how="inner", validate="one_to_one")
    merged[["ihc_status", "lost_markers"]] = merged.apply(
        lambda row: pd.Series(_classify_ihc(row)), axis=1
    )
    merged["msi_call"] = merged["msi_call"].str.upper()
    merged["msi_binary_call"] = merged["msi_call"].map({"MSI-H": 1, "MSS": 0}).fillna(0).astype(int)
    merged["ihc_dmmr"] = (merged["ihc_status"] == "dMMR").astype(int)

    merged["is_concordant"] = (
        ((merged["msi_call"] == "MSI-H") & (merged["ihc_status"] == "dMMR"))
        | ((merged["msi_call"] == "MSS") & (merged["ihc_status"] == "pMMR"))
    )
    merged["discordance_type"] = np.where(
        merged["is_concordant"],
        "Concordant",
        np.where(
            (merged["msi_call"] == "MSI-H") & (merged["ihc_status"] == "pMMR"),
            "MSI-H with retained IHC",
            "MSS with protein loss",
        ),
    )
    return merged


def summarize_concordance(merged_df: pd.DataFrame) -> Tuple[ConcordanceSummary, pd.DataFrame]:
    """Calculate concordance metrics."""

    total_cases = int(len(merged_df))
    concordant_cases = int(merged_df["is_concordant"].sum())
    discordant_cases = int(total_cases - concordant_cases)
    concordance_rate = concordant_cases / total_cases if total_cases else float("nan")
    discordance_rate = 1 - concordance_rate if total_cases else float("nan")

    contingency = pd.crosstab(merged_df["msi_call"], merged_df["ihc_status"])
    contingency_table = {
        str(msi_call): {str(ihc): int(count) for ihc, count in counts.items()}
        for msi_call, counts in contingency.iterrows()
    }

    discordant_df = merged_df.loc[~merged_df["is_concordant"], [
        "sample_id",
        "msi_call",
        "msi_score",
        "ihc_status",
        "lost_markers",
        "discordance_type",
    ]]

    summary = ConcordanceSummary(
        total_cases=total_cases,
        concordant_cases=concordant_cases,
        discordant_cases=discordant_cases,
        concordance_rate=concordance_rate,
        discordance_rate=discordance_rate,
        contingency_table=contingency_table,
        discordant_case_ids=discordant_df["sample_id"].tolist(),
    )
    return summary, discordant_df


def find_optimal_threshold(merged_df: pd.DataFrame) -> Dict[str, float]:
    """Identify the MSI score threshold that best matches IHC results."""

    scores = merged_df["msi_score"].to_numpy()
    labels = merged_df["ihc_dmmr"].to_numpy()
    thresholds = np.linspace(scores.min(), scores.max(), num=200)

    best = {
        "threshold": float("nan"),
        "youden_index": float("-inf"),
        "sensitivity": float("nan"),
        "specificity": float("nan"),
        "tp": 0,
        "tn": 0,
        "fp": 0,
        "fn": 0,
    }

    for threshold in thresholds:
        predicted = (scores >= threshold).astype(int)
        tp = int(((predicted == 1) & (labels == 1)).sum())
        tn = int(((predicted == 0) & (labels == 0)).sum())
        fp = int(((predicted == 1) & (labels == 0)).sum())
        fn = int(((predicted == 0) & (labels == 1)).sum())

        sensitivity = tp / (tp + fn) if (tp + fn) else float("nan")
        specificity = tn / (tn + fp) if (tn + fp) else float("nan")
        youden = (sensitivity + specificity - 1) if not np.isnan(sensitivity + specificity) else float("nan")

        if youden > best["youden_index"]:
            best.update(
                {
                    "threshold": float(threshold),
                    "youden_index": float(youden),
                    "sensitivity": float(sensitivity),
                    "specificity": float(specificity),
                    "tp": tp,
                    "tn": tn,
                    "fp": fp,
                    "fn": fn,
                }
            )

    return best


def fit_logistic_model(merged_df: pd.DataFrame) -> Tuple[Dict[str, float], pd.Series]:
    """Fit a logistic regression model linking MSI score to IHC dMMR status."""

    X = sm.add_constant(merged_df["msi_score"])
    y = merged_df["ihc_dmmr"]

    try:
        model = sm.Logit(y, X, missing="drop")
        result = model.fit(disp=False)
    except Exception as exc:  # pragma: no cover - defensive branch
        empty_predictions = pd.Series(np.nan, index=merged_df.index, name="predicted_prob_dmmr")
        return {"error": str(exc)}, empty_predictions

    conf_int = result.conf_int()
    logistic_summary = {
        "intercept": float(result.params["const"]),
        "msi_score_coef": float(result.params["msi_score"]),
        "p_values": {
            "intercept": float(result.pvalues["const"]),
            "msi_score": float(result.pvalues["msi_score"]),
        },
        "conf_int": {
            "intercept": [float(x) for x in conf_int.loc["const"].tolist()],
            "msi_score": [float(x) for x in conf_int.loc["msi_score"].tolist()],
        },
        "aic": float(result.aic),
    }

    predictions = pd.Series(result.predict(X), index=merged_df.index, name="predicted_prob_dmmr")
    return logistic_summary, predictions


def build_visualization(merged_df: pd.DataFrame, output_path: Path, threshold: float | None = None) -> None:
    """Create a plot summarising MSI score distribution by IHC status."""

    import matplotlib.pyplot as plt

    grouped = merged_df.groupby("ihc_status")
    labels = []
    data = []
    for status, subset in grouped:
        labels.append(status)
        data.append(subset["msi_score"].to_numpy())

    fig, ax = plt.subplots(figsize=(7, 5))
    box = ax.boxplot(data, labels=labels, patch_artist=True)
    colors = ["#c6d9f1", "#fde9d9"]
    for patch, color in zip(box["boxes"], colors[: len(box["boxes"])]):
        patch.set_facecolor(color)

    for idx, (status, subset) in enumerate(grouped, start=1):
        jitter = np.random.normal(loc=0, scale=0.04, size=len(subset))
        ax.scatter(
            np.full(len(subset), idx) + jitter,
            subset["msi_score"],
            alpha=0.7,
            edgecolor="black",
            linewidths=0.5,
            s=60,
            label=f"{status} samples",
        )

    if threshold is not None and np.isfinite(threshold):
        ax.axhline(threshold, color="#cc0000", linestyle="--", linewidth=1.5, label=f"Proposed threshold ({threshold:.1f})")

    ax.set_ylabel("MSI score")
    ax.set_title("MSI score distribution by IHC status")
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def run_analysis(msi_path: Path | str, ihc_path: Path | str, output_dir: Path | str) -> Dict[str, object]:
    """Execute the full concordance analysis pipeline."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    msi_df = load_msi_calls(msi_path)
    ihc_df = load_ihc_results(ihc_path)
    merged = merge_assay_results(msi_df, ihc_df)

    summary, discordant_df = summarize_concordance(merged)
    threshold_info = find_optimal_threshold(merged)
    logistic_info, predictions = fit_logistic_model(merged)

    merged = merged.copy()
    merged["predicted_prob_dmmr"] = predictions

    build_visualization(merged, output_dir / "msi_score_by_ihc_status.png", threshold_info.get("threshold"))

    merged.to_csv(output_dir / "merged_results.csv", index=False)
    discordant_df.to_csv(output_dir / "discordant_cases.csv", index=False)

    summary_dict = summary.to_dict()
    summary_dict["optimal_threshold"] = threshold_info
    summary_dict["logistic_regression"] = logistic_info

    summary_path = output_dir / "concordance_summary.json"
    summary_path.write_text(json.dumps(summary_dict, indent=2))

    return summary_dict

