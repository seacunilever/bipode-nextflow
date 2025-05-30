#!/usr/bin/env python3

import argparse
import logging
import os
from typing import Any, Dict, List, Optional, Tuple, Union

import multiqc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from multiqc.plots import scatter, table

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def calculate_summary_statistics(
    df: pd.Series,
) -> Dict[str, Union[np.ndarray, float, int]]:
    """Calculates summary statistics for BIFROST analysis."""
    logger.info("Calculating summary statistics...")

    # Cache frequently accessed values
    probes = np.array(df["probes"])
    max_conc = df["max_conc"]
    n_samp = df["n_samp"]
    conc = df["conc"]

    # Vectorize PoD mean calculation
    pod = np.array(
        [np.mean(df[i]["pod"]) if len(df[i]["pod"]) > 0 else np.nan for i in probes]
    )

    # Vectorize CDS calculation
    cds = np.array([df[i]["cds"] for i in probes])

    # Pre-allocate array and vectorize fold change calculation
    l2fc = np.empty(probes.shape[0], dtype="float")

    # Create a dictionary to cache response arrays
    response_cache = {probe: np.array(df[probe]["response"][1]) for probe in probes}

    # Vectorize fold change calculation
    for i, probe in enumerate(probes):
        y = response_cache[probe]
        index = np.argmax(np.abs(np.log2(y / y[0])))
        l2fc[i] = np.log2(y[index] / y[0])

    stats = {
        "probe": probes,
        "pod": pod,
        "cds": cds,
        "l2fc": l2fc,
        "max_conc": max_conc,
        "n_samp": n_samp,
        "conc": conc,
        "_response_cache": response_cache,
    }

    return stats


def filter_summary_statistics(
    df: Dict[str, np.ndarray], cds_threshold: float
) -> Dict[str, np.ndarray]:
    """Filters summary statistics based on CDS threshold.

    Args:
        df: Dictionary containing summary statistics with keys 'probe', 'pod', 'cds', 'l2fc'.
        cds_threshold: Minimum CDS value to keep in filtered results.

    Returns:
        Filtered dictionary with same structure as input, containing only entries where
        CDS >= cds_threshold.

    Raises:
        KeyError: If required keys are missing from input dictionary.
    """
    mask = df["cds"] >= cds_threshold
    filtered_df = df.copy()
    for key in ["probe", "pod", "cds", "l2fc"]:
        filtered_df[key] = df[key][mask]
    return filtered_df


def get_confidence_threshold_probability_density(
    x: np.ndarray,
    threshold_lower: float = 0.5,
    threshold_upper: float = 1.0,
    param_a: float = 0.38387606,
    param_b: float = -5.40387609,
    param_c: float = 2.8775016,
) -> np.ndarray:
    """Evaluates probability density for CDS threshold uncertainty.

    This function calculates the probability density for a defined function that describes
    uncertainty in the CDS threshold using a sigmoid-based transformation.

    Args:
        x: Array of values at which to calculate density.
        threshold_lower: Lower bound of the valid range (default: 0.5)
        threshold_upper: Upper bound of the valid range (default: 1.0)
        param_a: Parameter controlling the sigmoid shape (default: 0.38387606)
        param_b: Parameter controlling the sigmoid shift (default: -5.40387609)
        param_c: Parameter controlling the power transformation (default: 2.8775016)

    Returns:
        Array of probability density values corresponding to input x values.

    Note:
        The default parameters were determined empirically for CDS threshold uncertainty modeling.
        These values can be adjusted if different uncertainty modeling is required.
    """
    # Initialize output array with zeros
    probability_density = np.zeros(len(x))

    # Only calculate for values in the valid range
    valid_indices = np.where((x > threshold_lower) & (x < threshold_upper))[0]

    if len(valid_indices) == 0:
        return probability_density

    # Extract valid x values for cleaner calculations
    x_valid = x[valid_indices]

    # Step 1: Normalize x to [0,1] range and apply power transformation
    normalized_x = (x_valid - threshold_lower) / (threshold_upper - threshold_lower)
    g_function = normalized_x ** (-1 / param_c) - 1

    # Step 2: Calculate derivative of g with respect to x
    dg_dx = -(normalized_x ** (-1 / param_c - 1)) / (
        param_c * (threshold_upper - threshold_lower)
    )

    # Step 3: Apply logarithmic transformation
    h_function = param_b - np.log(g_function) / param_a

    # Step 4: Calculate derivative of h with respect to x
    dh_dx = -dg_dx / (param_a * g_function)

    # Step 5: Apply sigmoid function and its derivative to get probability density
    sigmoid_exp = np.exp(-h_function)
    sigmoid_derivative = sigmoid_exp / (1 + sigmoid_exp) ** 2
    probability_density[valid_indices] = sigmoid_derivative * dh_dx

    return probability_density


def get_minimum_pod_means(
    stats: Dict[str, np.ndarray], cds_thresholds: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes minimum PoD means for specified CDS thresholds.

    Args:
        stats: Dictionary containing probe statistics with keys 'probe', 'pod', 'cds', 'max_conc'.
        cds_thresholds: Array of CDS threshold values to evaluate.

    Returns:
        Tuple containing:
            - min_means: Array of minimum PoD means for each threshold
            - min_probes: Array of probe IDs corresponding to minimum means
            - min_cds: Array of CDS values for the minimum probes

    Raises:
        KeyError: If required keys are missing from stats dictionary.
    """
    min_means = np.full(len(cds_thresholds), stats["max_conc"])
    min_probes = np.full(len(cds_thresholds), "Max. conc.", dtype="object")
    min_cds = np.full(len(cds_thresholds), 0, dtype="float")

    for i, threshold in enumerate(cds_thresholds):
        mask = stats["cds"] >= threshold
        if np.sum(mask) > 0:
            pod = stats["pod"][mask]
            index = np.argmin(pod)
            min_means[i] = pod[index]
            min_probes[i] = stats["probe"][mask][index]
            min_cds[i] = stats["cds"][mask][index]

    return min_means, min_probes, min_cds


def get_global_pod(
    stats: Dict[str, np.ndarray],
) -> Dict[str, Union[float, np.ndarray, int]]:
    """Calculates global PoD from probe-level statistics.

    This function computes the global point of departure (PoD) by aggregating
    probe-level PoD distributions using a weighted approach based on CDS thresholds.

    Args:
        stats: Dictionary containing probe statistics with keys 'probe', 'pod', 'cds'.

    Returns:
        Dictionary containing:
            - global_pod: Calculated global PoD value
            - num_hits: Expected number of hits
            - means: Array of minimum PoD means for each threshold
            - probes: Array of probe IDs corresponding to minimum means
            - weights: Array of weights used in calculation
            - quantiles: Array of CDS threshold quantiles
            - cds: Array of CDS values for minimum probes

    Raises:
        KeyError: If required keys are missing from stats dictionary.
        ValueError: If stats data is malformed.
    """
    logger.info("Starting global PoD calculation...")

    # Get dictionary of minimum probes
    dq = 0.025
    quantiles = np.arange(0.5, 1 + dq, dq)
    quantiles[-1] -= 1e-6
    logger.info(f"Calculating minimum PoD means for {len(quantiles)} quantiles...")
    min_means, min_probes, min_cds = get_minimum_pod_means(stats, quantiles)

    # Calculate weights and global PoD
    logger.info("Calculating confidence threshold probability density...")
    weights = get_confidence_threshold_probability_density(quantiles)
    weight_sum = np.sum(weights)
    global_pod = 10 ** (np.sum(min_means * weights) / weight_sum)
    logger.info(f"Global PoD calculated: {global_pod}")

    # Calculate number of hits at each confidence threshold
    logger.info("Calculating number of hits at each confidence threshold...")
    num_hits = np.array([np.sum(stats["cds"] >= i) for i in quantiles])
    expected_num_hits = np.round(np.sum(num_hits * weights) / weight_sum)
    logger.info(f"Expected number of hits: {expected_num_hits}")

    results = {
        "global_pod": global_pod,
        "num_hits": expected_num_hits,
        "means": min_means,
        "probes": min_probes,
        "weights": weights,
        "quantiles": quantiles,
        "cds": min_cds,
    }

    return results


def get_min_probe_weights(df: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Computes weights for probes contributing to global PoD.

    Args:
        df: Dictionary containing global PoD data with keys 'probes', 'means', 'weights', 'cds'.

    Returns:
        Dictionary containing:
            - probe: Array of probe identifiers
            - weight: Array of weights for each probe
            - min_mean: Array of minimum means for each probe
            - cds: Array of CDS values for each probe

    Raises:
        KeyError: If required keys are missing from input dictionary.
    """
    # Aggregate weights for each unique probe
    probes = np.unique(df["probes"])
    means = np.power(10, [df["means"][df["probes"] == i][0] for i in probes])
    cds = np.array([df["cds"][df["probes"] == i][0] for i in probes])
    weights = np.array(
        [
            np.sum(df["weights"][df["probes"] == i]) / np.sum(df["weights"])
            for i in probes
        ]
    )

    rank = np.argsort(weights)[::-1]
    probes, means, cds, weights = probes[rank], means[rank], cds[rank], weights[rank]

    weight_dict = {"probe": probes, "weight": weights, "min_mean": means, "cds": cds}
    return weight_dict


def fit_pod_histogram(
    pod_samples: np.ndarray, n_samp: int
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Creates histogram approximation of PoD distribution.

    Args:
        pod_samples: Array of samples from the PoD distribution.
        n_samp: Maximum number of possible samples.

    Returns:
        Tuple containing:
            - weights: Array of histogram weights (or None if no valid samples)
            - bin_edges: Array of histogram bin edges (or None if no valid samples)

    Note:
        If pod_samples contains infinite values, n_samp is set to the length of pod_samples.
    """
    if np.isinf(pod_samples).any():
        n_samp = len(pod_samples)
    pod_samples = pod_samples[~np.isinf(pod_samples)]
    n_pod_samples = len(pod_samples)

    if n_pod_samples == 0:
        return None, None

    prob_response = n_pod_samples / n_samp
    counts, bin_edges = np.histogram(pod_samples, bins=int(np.sqrt(n_pod_samples)))
    bin_size = bin_edges[1:] - bin_edges[:-1]
    bin_weights = counts * bin_size
    total_weight = np.sum(bin_weights)
    weights = counts / total_weight * prob_response

    return weights, bin_edges


class ProbeData:
    """Helper class to manage probe data and calculations."""

    def __init__(self, df: pd.Series, probe: str, conc_units: str):
        self.df = df
        self.probe = probe
        self.conc_units = conc_units
        self._cache = {}
        # Pre-calculate commonly used values
        self._cache["cds"] = float(df[probe]["cds"])
        if self._cache["cds"] > 0:
            self._cache["mean_pod"] = np.mean(df[probe]["pod"])
            self._cache["pod_percentiles"] = self._calculate_pod_percentiles()

    def _calculate_pod_percentiles(
        self,
    ) -> Optional[Tuple[np.ndarray, List[int], List[float], List[str]]]:
        """Calculate PoD percentiles and related data if CDS > 0."""
        if self._cache["cds"] <= 0:
            return None
        percentiles = [1, 5, 10, 25, 75, 90, 95, 99]
        pod_percentiles = np.percentile(self.df[self.probe]["pod"], percentiles)
        pod_widths = [1, 1.5, 2, 2.5, 2.5, 2, 1.5, 1]
        pod_percentile_labels = [f"PoD {p}th percentile" for p in percentiles]
        return (pod_percentiles, percentiles, pod_widths, pod_percentile_labels)

    @property
    def cds(self) -> float:
        """Get CDS value for the probe."""
        return self._cache["cds"]

    @property
    def mean_pod(self) -> Optional[float]:
        """Calculate mean PoD if CDS > 0."""
        return self._cache.get("mean_pod")

    @property
    def pod_percentiles(
        self,
    ) -> Optional[Tuple[np.ndarray, List[int], List[float], List[str]]]:
        """Calculate PoD percentiles and related data if CDS > 0."""
        return self._cache.get("pod_percentiles")

    def get_response_data(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get response data for plotting."""
        if "response_data" not in self._cache:
            # Pre-calculate arrays to avoid repeated calculations
            conc = 10 ** np.array(
                self.df["conc"]
            )  # Convert log10 concentrations to actual values
            conc_index = np.array(self.df["conc_index"])
            count = np.array(self.df[self.probe]["count"])
            total_count = np.array(self.df["total_count"])
            median_total_count = np.median(total_count)

            # Use boolean masks for faster indexing
            treatment_mask = conc_index > 0
            control_mask = conc_index == 0

            # Calculate all values at once
            treatment_x = conc[conc_index[treatment_mask] - 1]
            treatment_y = (
                count[treatment_mask] / total_count[treatment_mask]
            ) * median_total_count
            control_y = (
                count[control_mask] / total_count[control_mask]
            ) * median_total_count
            response_x = 10 ** np.array(
                self.df[self.probe]["x"]
            )  # Convert log10 x values to actual values
            response = np.array(self.df[self.probe]["response"])

            self._cache["response_data"] = (
                treatment_x,
                treatment_y,
                control_y,
                response_x,
                response,
            )
        return self._cache["response_data"]


def create_summary_table_data(
    probes: List[str],
    df: pd.Series,
    stats: Dict[str, np.ndarray],
    weights: Dict[str, np.ndarray],
    conc_units: str,
    sort_by_abs_fc: bool = False,
) -> Dict[str, Dict[str, str]]:
    """Create summary table data for a list of probes.

    Args:
        probes: List of probe identifiers
        df: Pandas Series containing probe data
        stats: Dictionary containing summary statistics
        weights: Dictionary containing probe weights
        conc_units: String specifying concentration units
        sort_by_abs_fc: Whether to sort probes by absolute fold change

    Returns:
        Dictionary mapping probe IDs to their summary statistics
    """
    # First create all data without any sorting
    data = {}
    for probe in probes:
        # Calculate all probe statistics
        mean_pod = np.mean(df[probe]["pod"])
        weight = (
            weights["weight"][weights["probe"] == probe][0]
            if probe in weights["probe"]
            else 0.0
        )
        l2fc = stats["l2fc"][stats["probe"] == probe][0]

        # Store all data for this probe
        data[probe] = {
            "_abs_fc": abs(l2fc),  # Keep as float for potential sorting
            "CDS": f"{df[probe]['cds']:.3f}",
            "Mean PoD": f"{10**mean_pod:.2g} {conc_units}",
            "Log2 Fold Change": f"{l2fc:.2f}",
            "Global PoD Weight": f"{weight:.3f}",
            "Response Range": f"{df[probe]['response_threshold_lower']:.1f} - {df[probe]['response_threshold_upper']:.1f}",
        }

    # Then sort if requested, creating a new dictionary with sorted probes
    if sort_by_abs_fc:
        data = {
            probe: {**data[probe], "_abs_fc": f"{data[probe]['_abs_fc']:.3f}"}
            for probe in sorted(
                data.keys(), key=lambda x: data[x]["_abs_fc"], reverse=True
            )
        }
    else:
        # Just convert _abs_fc to string without sorting
        data = {
            probe: {**data[probe], "_abs_fc": f"{data[probe]['_abs_fc']:.3f}"}
            for probe in data
        }

    return data


def create_table_plot(
    data: Dict[str, Dict[str, str]],
    headers: Dict[str, Dict[str, Any]],
    table_id: str,
    title: str,
    sort_by_abs_fc: bool = False,
) -> table.plot:
    """Create a MultiQC table plot with common configuration."""
    pconfig = {
        "id": table_id,
        "title": title,
        "namespace": "BIFROST",
        "no_violin": True,
        "scale": False,  # Disable automatic scaling and coloring
        "sort_rows": False,  # Disable automatic sorting
        "col1_header": "Probe",  # This will label the first column as "Metric"
    }

    if sort_by_abs_fc:
        # Add _abs_fc to headers with supported options only
        headers["_abs_fc"] = {
            "title": "_abs_fc",
            "hidden": True,
            "description": "Absolute fold change (for sorting)",
            "placement": 0,  # Ensure it's the first column for sorting
        }

    return table.plot(data=data, headers=headers, pconfig=pconfig)


def fit_pod_histogram(
    pod_samples: np.ndarray, n_samp: int
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Creates histogram approximation of PoD distribution.

    Args:
        pod_samples: Array of samples from the PoD distribution.
        n_samp: Maximum number of possible samples.

    Returns:
        Tuple containing:
            - weights: Array of histogram weights (or None if no valid samples)
            - bin_edges: Array of histogram bin edges (or None if no valid samples)

    Note:
        If pod_samples contains infinite values, n_samp is set to the length of pod_samples.
    """
    if np.isinf(pod_samples).any():
        n_samp = len(pod_samples)
    pod_samples = pod_samples[~np.isinf(pod_samples)]
    n_pod_samples = len(pod_samples)

    if n_pod_samples == 0:
        return None, None

    prob_response = n_pod_samples / n_samp
    counts, bin_edges = np.histogram(pod_samples, bins=int(np.sqrt(n_pod_samples)))
    bin_size = bin_edges[1:] - bin_edges[:-1]
    bin_weights = counts * bin_size
    total_weight = np.sum(bin_weights)
    weights = counts / total_weight * prob_response

    return weights, bin_edges


def create_probe_plot(df: pd.Series, probe: str, conc_units: str) -> str:
    """Creates a concentration-response plot for a specific probe using Plotly."""
    logger.info(f"Creating concentration-response plot for probe {probe}")

    probe_data = ProbeData(df, probe, conc_units)
    treatment_x, treatment_y, control_y, response_x, response = (
        probe_data.get_response_data()
    )

    ymax = float(max(np.max(treatment_y), np.max(control_y), np.max(response[2])) * 1.1)

    fig = go.Figure()

    if df[probe]["cds"] > 0:
        mean_pod = probe_data.mean_pod
        weights, bin_edges = fit_pod_histogram(np.array(df[probe]["pod"]), df["n_samp"])
        if weights is not None and bin_edges is not None:
            bin_edges = 10**bin_edges
            weights = weights / np.max(weights) * 0.5 * df[probe]["cds"]

            shapes = []
            for i, weight in enumerate(weights):
                shapes.append(
                    dict(
                        type="rect",
                        line=dict(color="rgba(102, 51, 153, 0)"),
                        fillcolor=f"rgba(102, 51, 153, {weight})",
                        layer="below",
                        x0=bin_edges[i],
                        x1=bin_edges[i + 1],
                        xref="x",
                        y0=0,
                        y1=1,
                        yref="paper",
                    )
                )
            fig.update_layout(shapes=shapes)

    for i, y in enumerate(control_y):
        fig.add_hline(
            y=y,
            line=dict(color="#CCCCCC", width=0.6, dash="dash"),
            name="Solvent control" if i == 0 else None,
            showlegend=True if i == 0 else False,
            layer="below",
        )

    if df[probe]["cds"] > 0:
        fig.add_vline(
            x=float(10**mean_pod),
            line=dict(color="#663399", width=1.5, dash="solid"),
            name="Mean PoD | Response",
            showlegend=True,
        )

    fig.add_trace(
        go.Scatter(
            x=response_x,
            y=response[0],
            mode="lines",
            line=dict(color="#FF8080", width=1.5, dash="dash"),
            name="90% credible interval",
            showlegend=True,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=response_x,
            y=response[2],
            mode="lines",
            line=dict(color="#FF8080", width=1.5, dash="dash"),
            fill="tonexty",
            fillcolor="rgba(255, 128, 128, 0.1)",
            name=None,
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=response_x,
            y=response[1],
            mode="lines",
            line=dict(color="#FF0000", width=1.5),
            name="Median response",
            showlegend=True,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=treatment_x,
            y=treatment_y,
            mode="markers",
            marker=dict(symbol="x", size=7, color="#000000", line=dict(width=0.3)),
            name="Treatment data",
            showlegend=True,
            hovertemplate=(
                f"Concentration: %{{x:.2g}} {conc_units}<br>"
                + "Normalized count: %{y:.2f}<br>"
                + "<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title=dict(text=f"Concentration ({conc_units})", font=dict(size=12)),
            type="log",
            showgrid=True,
            gridcolor="#E5E5E5",
            gridwidth=1,
            showline=True,
            linewidth=1,
            linecolor="#000000",
            range=[
                np.log10(float(10 ** (df["conc"][0] - 1))),
                np.log10(float(10 ** (df["conc"][-1] + 1))),
            ],
        ),
        yaxis=dict(
            title=dict(text="Normalised count", font=dict(size=12)),
            showgrid=True,
            gridcolor="#E5E5E5",
            gridwidth=1,
            showline=True,
            linewidth=1,
            linecolor="#000000",
            range=[0, ymax],
        ),
        showlegend=True,
        legend=dict(
            yanchor="bottom",
            y=-0.3,
            xanchor="center",
            x=0.5,
            font=dict(size=11),
            bgcolor="rgba(255, 255, 255, 0.8)",
            orientation="h",
        ),
        margin=dict(l=60, r=20, t=40, b=80),
        height=400,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="Arial, sans-serif", size=11),
        uirevision=True,
        hovermode="closest",
        hoverdistance=10,
        spikedistance=10,
    )

    plot_html = fig.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={
            "responsive": True,
            "displayModeBar": True,
            "staticPlot": False,
            "showTips": True,
            "showLink": False,
            "displaylogo": False,
            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
            "toImageButtonOptions": {
                "format": "png",
                "filename": f"bifrost_plot_{probe}",
                "height": 400,
                "width": None,
                "scale": 2,
            },
        },
    )

    logger.info(f"Completed concentration-response plot for {probe}")
    return plot_html


def create_diagnostic_table_data(
    df: pd.Series,
    conc_units: str,
    cds_threshold: float,
    apply_cds_threshold: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """Create diagnostic table data for all probes."""
    diagnostic_data = {}
    for probe in df["probes"]:
        diag_text = df[probe]["diagnostics"]
        probe_data = ProbeData(df, probe, conc_units)

        # Skip probes that don't meet CDS threshold if filtering is enabled
        if apply_cds_threshold and probe_data.cds <= cds_threshold:
            continue

        # Parse individual checks
        checks = {
            "Treedepth": "✓" if "Treedepth satisfactory" in diag_text else "✗",
            "Divergences": "✓" if "No divergent transitions" in diag_text else "✗",
            "E-BFMI": "✓" if "E-BFMI satisfactory" in diag_text else "✗",
            "ESS": "✓" if "effective sample size satisfactory" in diag_text else "✗",
            "R-hat": (
                "✓"
                if "Rank-normalized split R-hat values satisfactory for all parameters"
                in diag_text
                else "✗"
            ),
        }

        # Calculate biological relevance score
        cds = probe_data.cds
        mean_pod = (
            probe_data.mean_pod if probe_data.mean_pod is not None else float("inf")
        )

        bio_score = 1 if cds > cds_threshold else 0
        if not np.isinf(mean_pod):
            bio_score += (df["max_conc"] - mean_pod) / df["max_conc"]

        # Extract R-hat parameters if present
        rhat_params = []
        if "R-hat greater than 1.01" in diag_text:
            start_idx = diag_text.find("greater than 1.01:") + len("greater than 1.01:")
            end_idx = diag_text.find("Such high values")
            if start_idx > 0 and end_idx > start_idx:
                params_text = diag_text[start_idx:end_idx].strip()
                rhat_params = [p.strip() for p in params_text.split() if p.strip()]

        # Check for regularization recommendation
        needs_regularization = (
            "You should consider regularizating your model with additional prior information or a more effective parameterization"
            in diag_text
        )

        # Format Mean PoD
        if not np.isnan(mean_pod):
            mean_pod_value = 10**mean_pod
            mean_pod_str = f"{mean_pod_value:.2g} {conc_units}"
        else:
            mean_pod_value = float("inf")
            mean_pod_str = "No response"

        # Add to diagnostic data
        diagnostic_data[probe] = {
            "CDS": float(cds),
            "CDS_str": f"{cds:.3f}",
            "Mean PoD": mean_pod_value,
            "Mean PoD_str": mean_pod_str,
            "Treedepth": checks["Treedepth"],
            "Divergences": checks["Divergences"],
            "E-BFMI": checks["E-BFMI"],
            "ESS": checks["ESS"],
            "R-hat": checks["R-hat"],
            "High R-hat Parameters": str(len(rhat_params)) if rhat_params else "0",
            "Response Range": f"{df[probe]['response_threshold_lower']:.1f} - {df[probe]['response_threshold_upper']:.1f}",
            "Needs Regularization": "⚠️" if needs_regularization else "✓",
            "_sort_score": bio_score,
        }

    return diagnostic_data


def get_plot_elements_description(
    cds_threshold: float, apply_cds_threshold: bool = False
) -> str:
    """Generate the common plot elements description used across multiple sections.

    Args:
        cds_threshold: The CDS threshold value
        apply_cds_threshold: Whether CDS threshold filtering is applied

    Returns:
        HTML string containing the plot elements description
    """
    return f"""
        <p><strong>Plot Elements:</strong></p>
        <ul>
            <li><strong>Black X markers</strong>: Normalised count for treated samples</li>
            <li><strong>Horizontal grey dashed lines</strong>: Normalised count for solvent control samples</li>
            <li><strong>Red bands</strong>: 90% credible interval for expected counts (light red fill)</li>
            <li><strong>Red line</strong>: Median of the expected count distribution</li>
            <li><strong>Purple bands</strong>: Histogram-like representation of PoD distribution (intensity indicates bin height)</li>
            <li><strong>Purple vertical line</strong>: Mean PoD estimate, conditional on there being a response</li>
            <li><strong>CDS</strong>: Concentration-dependency-score (probability of response below max concentration){" with threshold = " + str(cds_threshold) if apply_cds_threshold else ""}</li>
        </ul>
    """


def create_multiqc_report(
    summary_file,
    test_substance,
    cell_type,
    timepoint,
    conc_units,
    output_name,
    interactive_plots=False,
    n_fold_change_probes=5,
    cds_threshold=0.5,
    n_lowest_means=10,
    n_pod_stats=100,
    control_line_tolerance=0.02,
    min_control_lines=2,
    plot_height=400,
    pod_vs_fc_height=600,
    plots_force_flat_numseries=10000,
    no_cds_threshold=False,
):
    """Create a MultiQC report from BIFROST data."""
    logger.info(f"Starting report generation for {test_substance} on {cell_type}")

    # Convert no_cds_threshold to apply_cds_threshold (inverted logic)
    apply_cds_threshold = not no_cds_threshold

    # Configure MultiQC for interactive plots if requested
    if interactive_plots:
        os.environ["MULTIQC_PLOTS_FORCE_INTERACTIVE"] = "true"
        os.environ["MULTIQC_PLOTS_FLAT_NUMSERIES"] = str(plots_force_flat_numseries)

    # Load and process data
    logger.info(f"Loading summary data from {summary_file}")
    compression = "zip" if summary_file.endswith(".zip") else None
    df = pd.read_json(
        summary_file,
        typ="series",
        orient="index",
        compression=compression,
        dtype_backend="numpy_nullable",
    )

    # Compute summary statistics and global PoD
    stats = calculate_summary_statistics(df)
    stats_filtered = filter_summary_statistics(stats, cds_threshold=cds_threshold)
    global_pod = get_global_pod(stats)
    weights = get_min_probe_weights(global_pod)

    # Initialize MultiQC
    multiqc.reset()
    multiqc.config.plots_force_flat = not interactive_plots
    multiqc.config.plots_flat_numseries = (
        plots_force_flat_numseries if interactive_plots else 1000
    )
    multiqc.config.skip_generalstats = True
    multiqc.config.skip_plots = False
    multiqc.config.skip_cleanup = True

    # Verify interactive plots configuration after initialization
    if interactive_plots:
        multiqc.config.plots_force_interactive = True
        multiqc.config.plots_flat_numseries = plots_force_flat_numseries

    # Create summary table
    logger.info("Creating summary table...")
    # Create summary table - format data as a dictionary of samples (metrics)
    summary_data = {
        "Global PoD": {"Value": f"{global_pod['global_pod']:.2g} {conc_units}"},
        "Maximum tested concentration": {
            "Value": f"{10**stats['max_conc']:.2g} {conc_units}"
        },
        "Number of probes analyzed": {"Value": str(len(stats["probe"]))},
        "Number of hits": {"Value": str(int(global_pod["num_hits"]))},
        "Number of CDS>{cds_threshold} / CDS=1.0": {
            "Value": str(int(np.sum(stats["cds"] > cds_threshold)))
        },
        "Number of CDS=1.0": {"Value": str(int(np.sum(stats["cds"] == 1.0)))},
    }

    # Add minimum responding probe(s)
    valid_probes = weights["probe"][weights["probe"] != "Max. conc."]
    if len(valid_probes) > 0:
        # Sort by weight and get the top probe
        order = np.argsort(weights["weight"][weights["probe"] != "Max. conc."])[::-1]
        top_probe = valid_probes[order][0]
        top_weight = weights["weight"][weights["probe"] != "Max. conc."][order][0]
        top_cds = weights["cds"][weights["probe"] != "Max. conc."][order][0]
        top_pod = weights["min_mean"][weights["probe"] != "Max. conc."][order][0]

        summary_data["Minimum responding probe"] = {
            "Value": f"{top_probe}, weight={top_weight:.2g}, CDS={top_cds:.2g}, Mean PoD={top_pod:.2g} {conc_units}"
        }

    # Add largest fold changes
    if len(stats["l2fc"]) > 0:
        max_fc_idx = np.argmax(stats["l2fc"])
        min_fc_idx = np.argmin(stats["l2fc"])
        summary_data["Largest fold increase"] = {"Value": stats["probe"][max_fc_idx]}
        summary_data["Largest fold decrease"] = {"Value": stats["probe"][min_fc_idx]}

    # Add summary table to report
    summary_table = table.plot(
        data=summary_data,
        headers={"Value": {"title": "Value"}},
        pconfig={
            "id": "bifrost_summary",
            "title": "BIFROST Analysis Summary",
            "namespace": "BIFROST",
            "no_violin": True,
            "scale": False,  # Disable automatic scaling and coloring
            "sort_rows": False,
            "col1_header": "Metric",  # This will label the first column as "Metric"
        },
    )

    # Create PoD vs fold-change plot data
    pod_vs_fc_data = {
        "Probe": [
            {  # Empty string as dataset name to avoid showing in tooltips
                "x": float(x),  # Convert numpy float to Python float
                "y": float(y),  # Convert numpy float to Python float
                "text": str(text),  # Convert numpy string to Python string
                "name": str(text),  # Add name for hover text
            }
            for x, y, text in zip(
                10 ** stats_filtered["pod"],
                stats_filtered["l2fc"],
                stats_filtered["probe"],
            )
        ]
    }

    # Calculate y-axis limits
    if len(stats_filtered["l2fc"]) > 0:
        ymin, ymax = (
            min(stats_filtered["l2fc"].min(), -2) - 1,
            max(stats_filtered["l2fc"].max(), 2) + 1,
        )
    else:
        ymin, ymax = -2, 2

    # Calculate x-axis maximum
    xmax = (
        10 ** np.array(stats_filtered["conc"]).max() * 2
    )  # Double the max concentration for padding

    # Create scatter plot using MultiQC's scatter plot type
    pod_vs_fc_plot = scatter.plot(
        pod_vs_fc_data,
        pconfig={
            "id": "pod_vs_fc",
            "title": f'PoD vs Fold Change{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
            "xlab": f"Mean PoD | Response ({conc_units})",
            "ylab": "Max./min. log2 fold-change",
            "xlog": True,
            "xmin": 0,  # Set axis minimum to 0
            "xmax": xmax,  # Set axis maximum based on data
            "x_clipmin": 0.01,  # Clip data points below 0.01
            "x_clipmax": xmax,  # Clip data points above 100
            "x_decimals": 2,  # Format x-axis labels with 2 decimal places
            "ymin": ymin,  # Set y-axis minimum
            "ymax": ymax,  # Set y-axis maximum
            "marker_size": 5,
            "marker_line_width": 1,
            "color": "black",  # Use color instead of marker_line_color
            "opacity": 1.0,  # Set full opacity
            "showlegend": False,  # Hide legend
            "height": pod_vs_fc_height,  # Make plot taller to accommodate labels
            "x_lines": [  # Add vertical lines using x_lines
                {
                    "value": float(global_pod["global_pod"]),
                    "color": "#FF0000",  # Red
                    "width": 1,
                    "dash": "solid",
                    "label": "Global PoD",
                }
            ]
            + [
                {
                    "value": float(conc),
                    "color": "#D3D3D3",  # Light gray
                    "width": 1,
                    "dash": "dash",
                }
                for conc in 10 ** np.array(stats_filtered["conc"])
            ],
            "y_lines": [  # Add horizontal line at y=0
                {
                    "value": 0,
                    "color": "#CCCCCC",  # Light gray
                    "width": 1,
                    "dash": "dash",
                }
            ],
        },
    )

    # Create main BIFROST module for introduction and summary
    main_module = multiqc.BaseMultiqcModule(
        name="General",
        anchor="bifrost",
        href="https://github.com/your-repo/bifrost",
        info="BIFROST HTTr Analysis Report",
    )

    # Add introduction to main module
    main_module.add_section(
        name="Introduction",
        anchor="bifrost_intro",
        description=f"""
        <p>This report contains analysis of high-throughput transcriptomics data (HTTr) obtained after
        exposing {cell_type} cells for {timepoint} to {test_substance}. The BIFROST model
        (Bayesian inference for region of signal threshold) is a statistical model for analysis of HTTr concentration-response data.
        The model is designed to infer a point-of-departure (PoD) from a concentration-response dataset.
        The PoD is an estimate of the minimum effect concentration of the test substance
        for the experimental conditions under which the data were produced.
        PoDs are estimated as probability distributions.</p>

        <p>The implementation of the approach used here returns a single PoD for each probe analysed.
        PoD distributions are summarised in terms of quantiles of the distribution. The concentration-dependency-score (CDS) is the
        inferred probability that the test substance induces a change in expression below the maximum
        concentration tested.</p>

        <p>PoD distributions from individual probes are used to calculate a global PoD, defined as an estimate of a
        minimum effect concentration to induce perturbation in expression of any gene. The global PoD is formally
        an expectation with respect to the nominal concentration of the test substance.</p>

        <p>The report contains:</p>
        <ol>
            <li>A summary section including the global PoD and other overall statistics. Included within this section
            is a plot of the median PoD against the maximum log₂ fold-change in expression within the
            concentration-range.</li>
            <li>Concentration-response plots for probes with the 10 lowest expected PoDs. The entire probe set is
            first filtered for probes with a CDS > 0.5.</li>
            <li>A table summarising PoD statistics for the lowest 100 probes when ranked by the mean of the
            distribution conditional on there being a response.</li>
            <li>All PoDs are expressed with respect to the nominal concentration of the test substance.</li>
        </ol>
        """,
    )

    # Add summary section to main module
    main_module.add_section(
        name="Summary Statistics",
        anchor="bifrost_summary",
        plot=summary_table,
        description=f"""
        <p>Summary statistics from BIFROST analysis.</p>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li><strong>Global PoD</strong>: The estimated minimum effect concentration for any gene, summarizing the overall sensitivity of the system.</li>
            <li><strong>Maximum tested concentration</strong>: The highest concentration tested in the experiment.</li>
            <li><strong>Number of probes analyzed</strong>: Total number of probes included in the analysis.</li>
            <li><strong>Number of hits</strong>: Expected number of probes with a significant response.</li>
            <li><strong>Number of CDS>{cds_threshold} / CDS=1.0</strong>: Probes with strong concentration-dependent responses (CDS > {cds_threshold}) or maximal response (CDS = 1.0).</li>
            <li><strong>Minimum responding probe</strong>: The probe with the highest weight in the global PoD calculation, indicating the most sensitive response.</li>
            <li><strong>Largest fold increase/decrease</strong>: Probes with the largest positive or negative changes in expression.</li>
        </ul>
        """,
    )

    # Add PoD vs fold-change section to main module
    main_module.add_section(
        name="PoD vs Fold Change",
        anchor="bifrost_pod_vs_fc",
        plot=pod_vs_fc_plot,
        description=f"""
        <p>Maximum fold-change in expression over the tested concentration-range plotted against the
        probe-level PoD (mean given response){" for probes with CDS > " + str(cds_threshold) if apply_cds_threshold else ""}. The red vertical line indicates the global PoD.
        Vertical grey lines are placed at the experimental test substance concentrations.</p>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li>Each point represents a probe{" that meets the CDS threshold" if apply_cds_threshold else ""}.</li>
            <li><strong>X-axis (Mean PoD | Response)</strong>: Lower values indicate probes that respond at lower concentrations (more sensitive).</li>
            <li><strong>Y-axis (Max./min. log2 fold-change)</strong>: Higher absolute values indicate larger changes in expression.</li>
            <li>Probes in the lower-left region are most sensitive and show strong responses at low concentrations.</li>
            <li>The red vertical line (Global PoD) helps identify probes responding below the overall effect threshold.</li>
            {"<li>Only probes with CDS > " + str(cds_threshold) + " are shown, ensuring reliable concentration-dependent responses.</li>" if apply_cds_threshold else ""}
        </ul>
        """,
    )

    # Create module for weighted plots
    weighted_module = multiqc.BaseMultiqcModule(
        name=f'Probes with Non-zero Global PoD Weight{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
        anchor="bifrost_weighted",
        info=f'Concentration-response plots for probes contributing to global PoD{" (filtered by CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
    )

    # Create module for fold change plots
    fc_module = multiqc.BaseMultiqcModule(
        name=f'Fold Change Plots{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
        anchor="bifrost_fc_plots",
        info=f'PoD vs Fold Change plots{" for probes with CDS > " + str(cds_threshold) if apply_cds_threshold else ""}',
    )

    # Create module for lowest means plots
    lowest_means_module = multiqc.BaseMultiqcModule(
        name=f'Lowest Mean PoD Plots{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
        anchor="bifrost_lowest_means_plots",
        info=f'Concentration-response plots for probes with lowest mean PoD{" (filtered by CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
    )

    # Add description section to weighted module
    weighted_module.add_section(
        name="Overview",
        anchor="bifrost_weighted_overview",
        description=f"""
        <p>Concentration-response plots for probes that contribute to the global PoD calculation.</p>

        {get_plot_elements_description(cds_threshold, apply_cds_threshold)}

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>These probes were selected based on their contribution to the global PoD calculation.</li>
            <li>Each probe's weight in the global PoD calculation is determined by its CDS (Concentration-Dependency Score) and its position in the PoD distribution.</li>
            <li>Probes are sorted by their weight in descending order, showing the most influential probes first.</li>
            <li>Only probes with non-zero weights are included, as these are the ones that meaningfully contribute to the global PoD estimate.</li>
        </ul>
        """,
    )

    # Add plots for probes with non-zero global PoD weight to weighted module
    valid_probes = weights["probe"][weights["probe"] != "Max. conc."]
    probes_to_plot = valid_probes[
        np.argsort(weights["weight"][weights["probe"] != "Max. conc."])
    ]
    logger.info(f"Found {len(probes_to_plot)} probes with non-zero global PoD weight")

    if len(probes_to_plot) > 0:
        # Create summary table for weighted probes
        weighted_table_data = create_summary_table_data(
            probes_to_plot, df, stats, weights, conc_units, sort_by_abs_fc=True
        )
        weighted_summary_table = create_table_plot(
            data=weighted_table_data,
            headers={
                "CDS": {
                    "title": "CDS",
                    "description": "Concentration-Dependency Score",
                },
                "Mean PoD": {
                    "title": f"Mean PoD ({conc_units})",
                    "description": "Mean point of departure",
                },
                "Log2 Fold Change": {
                    "title": "Log2 Fold Change",
                    "description": "Maximum fold change in expression",
                },
                "Global PoD Weight": {
                    "title": "Global PoD Weight",
                    "description": "Weight in global PoD calculation",
                },
                "Response Range": {
                    "title": "Response Range",
                    "description": "Range of response thresholds",
                },
            },
            table_id="bifrost_weighted_summary",
            title="Summary Statistics for Probes with Non-zero Global PoD Weight",
            sort_by_abs_fc=True,
        )

        # Add summary table to weighted module
        weighted_module.add_section(
            name="Probe Summary Statistics",
            anchor="bifrost_weighted_summary",
            plot=weighted_summary_table,
            description="""
            <p>Summary statistics for probes contributing to the global PoD calculation.</p>
            <p>The table shows key metrics for each probe, sorted by their weight in the global PoD calculation.</p>
            """,
        )

        for probe in probes_to_plot:
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            weighted_module.add_section(
                name=probe,
                anchor=f"bifrost_weighted_{probe}",
                content=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}',
            )

    # Create module for probes with largest fold changes
    fc_module.add_section(
        name="Overview",
        anchor="bifrost_fc_overview",
        description=f"""
        <p>Concentration-response plots for probes with extreme expression changes{" that meet the CDS threshold" if apply_cds_threshold else ""}.</p>

        {get_plot_elements_description(cds_threshold, apply_cds_threshold)}

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>This section displays the {n_fold_change_probes} probes with the most extreme fold changes{" that meet two criteria:" if apply_cds_threshold else ":"}
                {"<ul>" if apply_cds_threshold else ""}
                {"<li>CDS > " + str(cds_threshold) + " (strong evidence for a concentration-dependent response)</li>" if apply_cds_threshold else ""}
                {"<li>Valid PoD estimate (mean PoD less than maximum tested concentration)</li>" if apply_cds_threshold else ""}
                {"</ul>" if apply_cds_threshold else ""}
                {"<li>Valid PoD estimate (mean PoD less than maximum tested concentration)</li>" if not apply_cds_threshold else ""}
            </li>
            {"<li>Probes are first filtered to include only those with CDS > " + str(cds_threshold) + ", ensuring reliable concentration-dependent responses.</li>" if apply_cds_threshold else ""}
            <li>Among {"these filtered probes" if apply_cds_threshold else "all probes"}, the {n_fold_change_probes} with the most extreme fold changes are selected.</li>
            <li>If fewer than {n_fold_change_probes} probes meet these criteria, all qualifying probes are shown.</li>
        </ul>
        """,
    )

    # Add section for most upregulated probes
    fc_module.add_section(
        name="Most Upregulated Probes",
        anchor="bifrost_fc_up",
        description=f"""
        <p>Concentration-response plots for the {n_fold_change_probes} probes with the largest positive fold changes (increased expression).</p>

        <p><strong>Selection Details:</strong></p>
        <ul>
            <li>These probes show the strongest increase in expression across the concentration range.</li>
            <li>Selected based on the maximum positive log2 fold change relative to control.</li>
            <li>Sorted by fold change magnitude in descending order.</li>
        </ul>
        """,
    )

    # Add plots for most upregulated probes
    if len(stats["l2fc"]) > 0:
        # Sort by absolute fold change magnitude
        abs_fc = np.abs(stats["l2fc"])
        index = np.argsort(abs_fc)[::-1]  # Sort in descending order
        n_up = min(n_fold_change_probes, len(stats["l2fc"]))
        # Get probes with largest absolute fold changes that are positive
        up_probes = stats["probe"][index][stats["l2fc"][index] > 0][:n_up]
        logger.info(f"Found {len(up_probes)} probes with largest positive fold changes")

        # Create summary table for upregulated probes
        up_table_data = create_summary_table_data(
            up_probes, df, stats, weights, conc_units, sort_by_abs_fc=True
        )
        up_summary_table = create_table_plot(
            data=up_table_data,
            headers={
                "CDS": {
                    "title": "CDS",
                    "description": "Concentration-Dependency Score",
                },
                "Mean PoD": {
                    "title": f"Mean PoD ({conc_units})",
                    "description": "Mean point of departure",
                },
                "Log2 Fold Change": {
                    "title": "Log2 Fold Change",
                    "description": "Maximum positive fold change in expression",
                },
                "Global PoD Weight": {
                    "title": "Global PoD Weight",
                    "description": "Weight in global PoD calculation",
                },
                "Response Range": {
                    "title": "Response Range",
                    "description": "Range of response thresholds",
                },
            },
            table_id="bifrost_fc_up_summary",
            title="Summary Statistics for Most Upregulated Probes",
            sort_by_abs_fc=True,
        )

        # Add summary table to upregulated section
        fc_module.add_section(
            name="Upregulated Probes Summary",
            anchor="bifrost_fc_up_summary",
            plot=up_summary_table,
            description="""
            <p>Summary statistics for the most upregulated probes.</p>
            <p>The table shows key metrics for each probe, sorted by their fold change magnitude.</p>
            """,
        )

        for probe in up_probes:
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            fc_module.add_section(
                name=probe,
                anchor=f"bifrost_fc_up_{probe}",
                content=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}, Log2 Fold Change = {stats["l2fc"][stats["probe"] == probe][0]:.2f}',
            )

    # Add section for most downregulated probes
    fc_module.add_section(
        name="Most Downregulated Probes",
        anchor="bifrost_fc_down",
        description=f"""
        <p>Concentration-response plots for the {n_fold_change_probes} probes with the largest negative fold changes (decreased expression).</p>

        <p><strong>Selection Details:</strong></p>
        <ul>
            <li>These probes show the strongest decrease in expression across the concentration range.</li>
            <li>Selected based on the maximum negative log2 fold change relative to control.</li>
            <li>Sorted by fold change magnitude in ascending order (most negative first).</li>
        </ul>
        """,
    )

    # Add plots for most downregulated probes
    if len(stats["l2fc"]) > 0:
        # Sort by absolute fold change magnitude
        abs_fc = np.abs(stats["l2fc"])
        index = np.argsort(abs_fc)[::-1]  # Sort in descending order
        n_down = min(n_fold_change_probes, len(stats["l2fc"]))
        # Get probes with largest absolute fold changes that are negative
        down_probes = stats["probe"][index][stats["l2fc"][index] < 0][:n_down]
        logger.info(f"Found {len(down_probes)} probes with largest negative fold changes")

        # Create summary table for downregulated probes
        down_table_data = create_summary_table_data(
            down_probes, df, stats, weights, conc_units, sort_by_abs_fc=True
        )
        down_summary_table = create_table_plot(
            data=down_table_data,
            headers={
                "CDS": {
                    "title": "CDS",
                    "description": "Concentration-Dependency Score",
                },
                "Mean PoD": {
                    "title": f"Mean PoD ({conc_units})",
                    "description": "Mean point of departure",
                },
                "Log2 Fold Change": {
                    "title": "Log2 Fold Change",
                    "description": "Maximum negative fold change in expression",
                },
                "Global PoD Weight": {
                    "title": "Global PoD Weight",
                    "description": "Weight in global PoD calculation",
                },
                "Response Range": {
                    "title": "Response Range",
                    "description": "Range of response thresholds",
                },
            },
            table_id="bifrost_fc_down_summary",
            title="Summary Statistics for Most Downregulated Probes",
            sort_by_abs_fc=True,
        )

        # Add summary table to downregulated section
        fc_module.add_section(
            name="Downregulated Probes Summary",
            anchor="bifrost_fc_down_summary",
            plot=down_summary_table,
            description="""
            <p>Summary statistics for the most downregulated probes.</p>
            <p>The table shows key metrics for each probe, sorted by their fold change magnitude.</p>
            """,
        )

        for probe in down_probes:
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            fc_module.add_section(
                name=probe,
                anchor=f"bifrost_fc_down_{probe}",
                content=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}, Log2 Fold Change = {stats["l2fc"][stats["probe"] == probe][0]:.2f}',
            )

    # Add plots for lowest means to lowest means module
    n_probe = len(stats["probe"])
    if apply_cds_threshold:
        mask = stats["cds"] > cds_threshold
        probes_to_plot = stats["probe"][mask][np.argsort(stats["pod"][mask])][
            : min(np.sum(mask), n_lowest_means)
        ]
    else:
        probes_to_plot = stats["probe"][np.argsort(stats["pod"])][
            : min(n_probe, n_lowest_means)
        ]
    logger.info(f"Found {len(probes_to_plot)} probes with lowest means to plot")

    if len(probes_to_plot) > 0:
        # Add overview section to lowest means module
        lowest_means_module.add_section(
            name="Overview",
            anchor="bifrost_lowest_means_overview",
            description=f"""
        {get_plot_elements_description(cds_threshold, apply_cds_threshold)}

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>This section displays the {n_lowest_means} probes with the lowest mean PoD values{" that meet two criteria:" if apply_cds_threshold else ":"}
                {"<ul>" if apply_cds_threshold else ""}
                {"<li>CDS > " + str(cds_threshold) + " (strong evidence for a concentration-dependent response)</li>" if apply_cds_threshold else ""}
                {"<li>Valid PoD estimate (mean PoD less than maximum tested concentration)</li>" if apply_cds_threshold else ""}
                {"</ul>" if apply_cds_threshold else ""}
                {"<li>Valid PoD estimate (mean PoD less than maximum tested concentration)</li>" if not apply_cds_threshold else ""}
            </li>
            {"<li>Probes are first filtered to include only those with CDS > " + str(cds_threshold) + ", ensuring reliable concentration-dependent responses.</li>" if apply_cds_threshold else ""}
            <li>Among {"these filtered probes" if apply_cds_threshold else "all probes"}, the {n_lowest_means} with the lowest mean PoD values are selected.</li>
            <li>If fewer than {n_lowest_means} probes meet these criteria, all qualifying probes are shown.</li>
        </ul>
        """,
        )

        # Create summary table for lowest means probes
        lowest_means_table_data = create_summary_table_data(
            probes_to_plot, df, stats, weights, conc_units, sort_by_abs_fc=True
        )
        lowest_means_summary_table = create_table_plot(
            data=lowest_means_table_data,
            headers={
                "CDS": {
                    "title": "CDS",
                    "description": "Concentration-Dependency Score",
                },
                "Mean PoD": {
                    "title": f"Mean PoD ({conc_units})",
                    "description": "Mean point of departure",
                },
                "Log2 Fold Change": {
                    "title": "Log2 Fold Change",
                    "description": "Maximum fold change in expression",
                },
                "Global PoD Weight": {
                    "title": "Global PoD Weight",
                    "description": "Weight in global PoD calculation",
                },
                "Response Range": {
                    "title": "Response Range",
                    "description": "Range of response thresholds",
                },
            },
            table_id="bifrost_lowest_means_summary",
            title="Summary Statistics for Most Sensitive Probes (CDS > 0.5)",
            sort_by_abs_fc=True,
        )

        # Add summary table to lowest means module
        lowest_means_module.add_section(
            name="Probe Summary Statistics",
            anchor="bifrost_lowest_means_summary",
            plot=lowest_means_summary_table,
            description=f"""
            <p>Summary statistics for the most sensitive probes with CDS > 0.5.</p>
            <p>The table shows key metrics for each probe, sorted by their mean PoD (most sensitive first).</p>
            """,
        )

        for probe in probes_to_plot:
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            lowest_means_module.add_section(
                name=probe,
                anchor=f"bifrost_lowest_means_{probe}",
                content=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}',
            )

    # Create module for PoD statistics
    stats_module = multiqc.BaseMultiqcModule(
        name=f'Probe-level PoD Statistics{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
        anchor="bifrost_stats",
        info=f'Detailed statistics for probes{" with CDS > " + str(cds_threshold) if apply_cds_threshold else ""}',
    )

    # Add PoD statistics table to stats module
    if n_probe > 0:
        # Get probes and sort by PoD
        if apply_cds_threshold:
            cds_mask = stats["cds"] > cds_threshold
            probes = stats["probe"][cds_mask][np.argsort(stats["pod"][cds_mask])][
                :n_pod_stats
            ]
            cds = stats["cds"][cds_mask][np.argsort(stats["pod"][cds_mask])][
                :n_pod_stats
            ]
        else:
            probes = stats["probe"][np.argsort(stats["pod"])][:n_pod_stats]
            cds = stats["cds"][np.argsort(stats["pod"])][:n_pod_stats]
        logger.info(f"Found {len(probes)} probes to include in PoD statistics table")

        # Create table data
        table_data = {}
        for probe, cds_val in zip(probes, cds):
            n_pod_samples = len(df[probe]["pod"])
            extended_pod_samples = np.concatenate(
                (
                    df[probe]["pod"],
                    [df["max_conc"] for _ in range(df["n_samp"] - n_pod_samples)],
                )
            )
            pod_percentiles = np.percentile(extended_pod_samples, q=(5, 25, 50, 75, 95))

            # Format PoD values
            pod_values = []
            for pod_val in pod_percentiles:
                if pod_val < df["max_conc"]:
                    # Convert to integer and format without scientific notation
                    pod_values.append(f"{int(10**pod_val)}")
                else:
                    pod_values.append(f">{int(10**pod_val)}")

            # Add row to table data
            table_data[probe] = {
                "CDS": f"{cds_val:.3f}",
                "5th percentile": pod_values[0],
                "25th percentile": pod_values[1],
                "50th percentile": pod_values[2],
                "75th percentile": pod_values[3],
                "95th percentile": pod_values[4],
            }

        # Create table plot
        pod_stats_table = table.plot(
            data=table_data,
            headers={
                "CDS": {
                    "title": "CDS",
                    "format": "{:.3f}",
                    "description": f'Concentration-Dependency Score{" (threshold = " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
                },
                "5th percentile": {
                    "title": f"5th percentile ({conc_units})",
                    "description": "5th percentile of PoD distribution",
                },
                "25th percentile": {
                    "title": f"25th percentile ({conc_units})",
                    "description": "25th percentile of PoD distribution",
                },
                "50th percentile": {
                    "title": f"50th percentile ({conc_units})",
                    "description": "Median of PoD distribution",
                },
                "75th percentile": {
                    "title": f"75th percentile ({conc_units})",
                    "description": "75th percentile of PoD distribution",
                },
                "95th percentile": {
                    "title": f"95th percentile ({conc_units})",
                    "description": "95th percentile of PoD distribution",
                },
            },
            pconfig={
                "id": "bifrost_stats_table",
                "title": f'Probe-level PoD Statistics{" (CDS > " + str(cds_threshold) + ")" if apply_cds_threshold else ""}',
                "namespace": "BIFROST",
                "no_violin": True,
                "scale": False,  # Disable automatic scaling and coloring
                "sort_rows": False,
                "col1_header": "Probe",  # Label first column as Probe
            },
        )

        # Add table to section
        stats_module.add_section(
            name="PoD Statistics Table",
            anchor="bifrost_stats_table",
            plot=pod_stats_table,
            description=f"""
            <p>Summary statistics for probes{" with CDS > " + str(cds_threshold) if apply_cds_threshold else ""}, showing the distribution of PoD (Point of Departure) values.</p>
            <p>The table includes:</p>
            <ul>
                <li><strong>CDS</strong>: Concentration-Dependency Score (probability of response below max concentration){" with threshold = " + str(cds_threshold) if apply_cds_threshold else ""}</li>
                <li><strong>PoD percentiles</strong>: Different quantiles of the PoD distribution</li>
                <li>Values are shown in {conc_units}</li>
                <li>Probes are sorted by median PoD (50th percentile)</li>
                <li>Only shows the top {n_pod_stats} probes{" with CDS > " + str(cds_threshold) if apply_cds_threshold else ""}</li>
            </ul>

            <p>Hover over column headers for more detailed descriptions.</p>
            """,
        )

    # Create diagnostic table data with parsed checks
    diagnostic_data = create_diagnostic_table_data(
        df, conc_units, cds_threshold, apply_cds_threshold
    )

    # Create diagnostic table
    diagnostic_table = create_table_plot(
        data={
            k: {
                sk: v[sk]
                for sk in [
                    "CDS_str",
                    "Mean PoD_str",
                    "Treedepth",
                    "Divergences",
                    "E-BFMI",
                    "ESS",
                    "R-hat",
                    "High R-hat Parameters",
                    "Response Range",
                    "Needs Regularization",
                    "_sort_score",
                ]
            }
            for k, v in diagnostic_data.items()
        },
        headers={
            "CDS_str": {
                "title": "CDS",
                "format": "{:.3f}",
                "description": f"Concentration-Dependency Score (probability of response below max concentration, threshold = {cds_threshold})",
                "cond_formatting_rules": {
                    "pass": [
                        {"gt": cds_threshold}
                    ]  # Highlight probes with CDS > threshold
                },
            },
            "Mean PoD_str": {
                "title": f"Mean PoD ({conc_units})",
                "description": 'Mean point of departure (effect concentration). "No response" indicates no valid PoD samples.',
                "cond_formatting_rules": {
                    "warn": [
                        {"s_eq": "No response"}
                    ]  # Highlight probes with no response
                },
            },
            "Treedepth": {
                "title": "Treedepth",
                "description": "Sampler transitions treedepth check",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for pass
                    "fail": [{"s_eq": "✗"}],  # Red for fail
                },
            },
            "Divergences": {
                "title": "Divergences",
                "description": "Check for divergent transitions",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for pass
                    "fail": [{"s_eq": "✗"}],  # Red for fail
                },
            },
            "E-BFMI": {
                "title": "E-BFMI",
                "description": "HMC potential energy check",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for pass
                    "fail": [{"s_eq": "✗"}],  # Red for fail
                },
            },
            "ESS": {
                "title": "ESS",
                "description": "Effective sample size check",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for pass
                    "fail": [{"s_eq": "✗"}],  # Red for fail
                },
            },
            "R-hat": {
                "title": "R-hat",
                "description": "Gelman-Rubin convergence diagnostic",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for pass
                    "fail": [{"s_eq": "✗"}],  # Red for fail
                },
            },
            "High R-hat Parameters": {
                "title": "# Parameters with R-hat > 1.01",
                "description": "Number of parameters with high R-hat values",
            },
            "Response Range": {
                "title": "Response Range",
                "description": "Range of response thresholds",
            },
            "Needs Regularization": {
                "title": "⚠️ Regularization",
                "description": "Model may need regularization",
                "cond_formatting_rules": {
                    "pass": [{"s_eq": "✓"}],  # Green for no regularization needed
                    "warn": [{"s_eq": "⚠️"}],  # Orange for needs regularization
                },
            },
            "_sort_score": {
                "title": "_sort_score",
                "hidden": True,  # Hide the sorting column from display
            },
        },
        table_id="bifrost_diagnostics_table",
        title="Probe Diagnostic Summary",
    )

    # Create module for diagnostics
    diag_module = multiqc.BaseMultiqcModule(
        name="Diagnostic Summary",
        anchor="bifrost_diagnostics",
        info="Model diagnostics and quality checks",
    )

    # Add diagnostic table to diag module
    diag_module.add_section(
        name="Diagnostic Table",
        anchor="bifrost_diagnostics_table",
        plot=diagnostic_table,
        description=f"""
        <p>Summary of diagnostic checks for each probe. The table is sorted to show the most biologically relevant probes first, based on:</p>
        <ul>
            <li><strong>CDS > {cds_threshold}</strong>: Probes with strong concentration-dependent responses are prioritized</li>
            <li><strong>Lower Mean PoD</strong>: Among probes with CDS > {cds_threshold}, those with lower PoD (more sensitive) are shown first</li>
        </ul>

        <p><strong>Model Convergence Checks</strong> (✓ = pass, ✗ = fail):</p>
        <ul>
            <li><strong>Treedepth</strong>: Checks if sampler transitions reached maximum treedepth</li>
            <li><strong>Divergences</strong>: Checks for divergent transitions in the sampler</li>
            <li><strong>E-BFMI</strong>: Checks Hamiltonian Monte Carlo potential energy</li>
            <li><strong>ESS</strong>: Checks effective sample size for all parameters</li>
            <li><strong>R-hat</strong>: Checks if all parameters have satisfactory rank-normalized split R-hat values</li>
        </ul>

        <p><strong>Additional Information:</strong></p>
        <ul>
            <li><strong>High R-hat Parameters</strong>: Number of parameters with R-hat > 1.01 that may indicate incomplete mixing</li>
            <li><strong>CDS</strong>: Concentration-Dependency Score (probability of response below max concentration)</li>
            <li><strong>Mean PoD</strong>: Mean point of departure (effect concentration)</li>
            <li><strong>Response Range</strong>: Range of response thresholds</li>
            <li><strong>Regularization</strong>: ⚠️ indicates model may need regularization with additional prior information or more effective parameterization</li>
        </ul>

        <p>Probes with "No response" in the Mean PoD column did not show a significant response in the tested range.
        Failed diagnostic checks (red ✗) indicate potential model issues for that probe. High R-hat values (>1.01) suggest
        that the model may need additional regularization or reparameterization to improve mixing.</p>
        """,
    )

    # Add all modules to report
    multiqc.report.modules.extend(
        [
            main_module,
            weighted_module,
            fc_module,
            lowest_means_module,
            stats_module,
            diag_module,
        ]
    )

    # Write report
    logger.info("Generating report...")

    try:
        multiqc.config.verbose = True
        multiqc.write_report(
            output_dir=os.path.dirname(output_name),
            filename=os.path.basename(output_name),
            title=f"BIFROST HTTr Analysis - {test_substance} ({cell_type})",
            report_comment=f"Analysis of {test_substance} on {cell_type} cells after {timepoint} exposure",
            force=True,
        )
    except Exception as e:
        logger.error(f"Error during report generation: {str(e)}")
        raise

    logger.info("Report generation complete")
    if not interactive_plots:
        logger.info(
            "Note: Consider using --interactive-plots for faster rendering with large datasets"
        )


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create Bifrost HTTR reports using MultiQC"
    )
    parser.add_argument(
        "--summary-file", required=True, help="Path to the summary JSON file"
    )
    parser.add_argument(
        "--test-substance",
        default="MyChemical",
        help="Name of the test substance (default: MyChemical)",
    )
    parser.add_argument(
        "--cell-type",
        default="MyCell",
        help="Type of cell used in the test (default: MyCell)",
    )
    parser.add_argument(
        "--output-name",
        default="multiqc_report.html",
        help="Name for the output report (default: multiqc_report.html)",
    )
    parser.add_argument(
        "--timepoint",
        default="24 hours",
        help="Exposure duration within experiment (default: 24 hours)",
    )
    parser.add_argument(
        "--conc-units",
        default="uM",
        choices=["uM", "ugml-1", "mgml-1"],
        help="Concentration units (default: uM)",
    )
    parser.add_argument(
        "--interactive-plots",
        action="store_true",
        help="Force interactive plots (may be faster for large datasets)",
    )
    parser.add_argument(
        "--n-fold-change-probes",
        type=int,
        default=5,
        help="Number of most up/down regulated probes to show (default: 5)",
    )
    parser.add_argument(
        "--cds-threshold",
        type=float,
        default=0.5,
        help="Concentration-Dependency Score threshold for filtering probes (default: 0.5)",
    )
    parser.add_argument(
        "--n-lowest-means",
        type=int,
        default=10,
        help="Number of lowest mean PoD probes to show (default: 10)",
    )
    parser.add_argument(
        "--n-pod-stats",
        type=int,
        default=100,
        help="Number of probes to include in PoD statistics table (default: 100)",
    )
    parser.add_argument(
        "--control-line-tolerance",
        type=float,
        default=0.02,
        help="Tolerance for filtering similar control lines (default: 0.02)",
    )
    parser.add_argument(
        "--min-control-lines",
        type=int,
        default=2,
        help="Minimum number of control lines to show (default: 2)",
    )
    parser.add_argument(
        "--plot-height",
        type=int,
        default=400,
        help="Height of concentration-response plots in pixels (default: 400)",
    )
    parser.add_argument(
        "--pod-vs-fc-height",
        type=int,
        default=600,
        help="Height of PoD vs Fold Change plot in pixels (default: 600)",
    )
    parser.add_argument(
        "--plots-force-flat-numseries",
        type=int,
        default=10000,
        help="Maximum number of series for flat plots (default: 10000)",
    )
    parser.add_argument(
        "--no-cds-threshold",
        action="store_true",
        help="Do not filter probes by CDS threshold in summary tables and lowest mean PoDs section (default: False, meaning filtering is applied)",
    )

    return parser.parse_args()


def main():
    """Main entry point for report generation script."""
    args = parse_args()

    create_multiqc_report(
        args.summary_file,
        args.test_substance,
        args.cell_type,
        args.timepoint,
        args.conc_units,
        args.output_name,
        args.interactive_plots,
        args.n_fold_change_probes,
        args.cds_threshold,
        args.n_lowest_means,
        args.n_pod_stats,
        args.control_line_tolerance,
        args.min_control_lines,
        args.plot_height,
        args.pod_vs_fc_height,
        args.plots_force_flat_numseries,
        args.no_cds_threshold,
    )


if __name__ == "__main__":
    main()
