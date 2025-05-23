#!/usr/bin/env python3

import os
import logging
import argparse
import json
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import multiqc
from multiqc.plots import table, linegraph, scatter, box
import time
import signal
from contextlib import contextmanager
import threading
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def calculate_summary_statistics(df: pd.Series) -> Dict[str, Union[np.ndarray, float, int]]:
    """Calculates summary statistics for BIFROST analysis."""
    logger.info("Calculating summary statistics...")

    # Cache frequently accessed values
    probes = np.array(df['probes'])
    max_conc = df['max_conc']
    n_samp = df['n_samp']
    conc = df['conc']

    # Vectorize PoD mean calculation
    pod = np.array([np.mean(df[i]['pod']) if len(df[i]['pod']) > 0 else np.nan for i in probes])

    # Vectorize CDS calculation
    cds = np.array([df[i]['cds'] for i in probes])

    # Pre-allocate array and vectorize fold change calculation
    l2fc = np.empty(probes.shape[0], dtype='float')

    # Create a dictionary to cache response arrays
    response_cache = {probe: np.array(df[probe]['response'][1]) for probe in probes}

    # Vectorize fold change calculation
    for i, probe in enumerate(probes):
        y = response_cache[probe]
        index = np.argmax(np.abs(np.log2(y / y[0])))
        l2fc[i] = np.log2(y[index] / y[0])

    stats = {
        'probe': probes,
        'pod': pod,
        'cds': cds,
        'l2fc': l2fc,
        'max_conc': max_conc,
        'n_samp': n_samp,
        'conc': conc,
        '_response_cache': response_cache
    }

    return stats

def filter_summary_statistics(df: Dict[str, np.ndarray], cds_threshold: float) -> Dict[str, np.ndarray]:
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
    mask = df['cds'] >= cds_threshold
    filtered_df = df.copy()
    for key in ['probe', 'pod', 'cds', 'l2fc']:
        filtered_df[key] = df[key][mask]
    return filtered_df

def get_confidence_threshold_probability_density(x: np.ndarray) -> np.ndarray:
    """Evaluates probability density for CDS threshold uncertainty.

    This function calculates the probability density for a defined function that describes
    uncertainty in the CDS threshold.

    Args:
        x: Array of values at which to calculate density.

    Returns:
        Array of probability density values corresponding to input x values.

    Note:
        The function uses predefined parameters (tl, tu, a, b, c) that were determined
        empirically for CDS threshold uncertainty modeling.
    """
    dq = np.zeros(len(x))
    index = np.where((x > 0.5) & (x < 1))[0]
    tl, tu, a, b, c = 0.5, 1, 0.38387606, -5.40387609, 2.8775016
    g = ((x[index] - tl) / (tu - tl)) ** (-1 / c) - 1
    dg = - ((x[index] - tl) / (tu - tl)) ** (-1 / c - 1) / (c * (tu - tl))
    h = b - np.log(g) / a
    dh = - dg / (a * g)
    dq[index] = np.exp(-h) / (1 + np.exp(-h)) ** 2 * dh
    return dq

def get_minimum_pod_means(stats: Dict[str, np.ndarray], cds_thresholds: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    min_means = np.full(len(cds_thresholds), stats['max_conc'])
    min_probes = np.full(len(cds_thresholds), 'Max. conc.', dtype='object')
    min_cds = np.full(len(cds_thresholds), 0, dtype='float')

    for i, threshold in enumerate(cds_thresholds):
        mask = stats['cds'] >= threshold
        if np.sum(mask) > 0:
            pod = stats['pod'][mask]
            index = np.argmin(pod)
            min_means[i] = pod[index]
            min_probes[i] = stats['probe'][mask][index]
            min_cds[i] = stats['cds'][mask][index]

    return min_means, min_probes, min_cds

def get_global_pod(stats: Dict[str, np.ndarray]) -> Dict[str, Union[float, np.ndarray, int]]:
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
    num_hits = np.array([np.sum(stats['cds'] >= i) for i in quantiles])
    expected_num_hits = np.round(np.sum(num_hits * weights) / weight_sum)
    logger.info(f"Expected number of hits: {expected_num_hits}")

    results = {
        'global_pod': global_pod,
        'num_hits': expected_num_hits,
        'means': min_means,
        'probes': min_probes,
        'weights': weights,
        'quantiles': quantiles,
        'cds': min_cds
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
    probes = np.unique(df['probes'])
    means = np.power(10, [df['means'][df['probes'] == i][0] for i in probes])
    cds = np.array([df['cds'][df['probes'] == i][0] for i in probes])
    weights = np.array([np.sum(df['weights'][df['probes'] == i]) / np.sum(df['weights']) for i in probes])

    rank = np.argsort(weights)[::-1]
    probes, means, cds, weights = probes[rank], means[rank], cds[rank], weights[rank]

    weight_dict = {
        'probe': probes,
        'weight': weights,
        'min_mean': means,
        'cds': cds
    }
    return weight_dict

def filter_similar_control_lines(control_y: np.ndarray, tolerance: float = 0.02, min_lines: int = 2) -> np.ndarray:
    """Filter out control lines that are too similar to each other.

    Args:
        control_y: Array of control y-values
        tolerance: Maximum relative difference between y-values to consider them similar (default: 0.02 or 2%)
        min_lines: Minimum number of control lines to show (default: 2)

    Returns:
        Filtered array of control y-values
    """
    if len(control_y) <= min_lines:
        return control_y

    # Sort y values
    sorted_y = np.sort(control_y)

    # Keep first value
    filtered = [sorted_y[0]]

    # Check each subsequent value against the last kept value
    for y in sorted_y[1:]:
        # Handle zero values
        if filtered[-1] == 0:
            if y == 0:
                continue  # Skip if both values are zero
            else:
                filtered.append(y)  # Keep non-zero value
        else:
            # Calculate relative difference for non-zero values
            rel_diff = abs(y - filtered[-1]) / filtered[-1]
            if rel_diff > tolerance:
                filtered.append(y)

    # If we have fewer than min_lines, add more values
    if len(filtered) < min_lines:
        # Get remaining values that weren't included
        remaining = sorted_y[~np.isin(sorted_y, filtered)]
        # Add values until we reach min_lines or run out of values
        for y in remaining:
            if len(filtered) >= min_lines:
                break
            filtered.append(y)
        # Sort the final list to maintain order
        filtered.sort()

    return np.array(filtered)

class ProbeData:
    """Helper class to manage probe data and calculations."""
    def __init__(self, df: pd.Series, probe: str, conc_units: str):
        self.df = df
        self.probe = probe
        self.conc_units = conc_units
        self._cache = {}

    @property
    def cds(self) -> float:
        """Get CDS value for the probe."""
        return float(self.df[self.probe]['cds'])

    @property
    def mean_pod(self) -> Optional[float]:
        """Calculate mean PoD if CDS > 0."""
        if self.cds <= 0:
            return None
        if 'mean_pod' not in self._cache:
            self._cache['mean_pod'] = np.mean(self.df[self.probe]['pod'])
        return self._cache['mean_pod']

    @property
    def pod_percentiles(self) -> Optional[Tuple[np.ndarray, List[int], List[float], List[str]]]:
        """Calculate PoD percentiles and related data if CDS > 0."""
        if self.cds <= 0:
            return None
        if 'pod_percentiles' not in self._cache:
            percentiles = [1, 5, 10, 25, 75, 90, 95, 99]
            pod_percentiles = np.percentile(self.df[self.probe]['pod'], percentiles)
            pod_widths = [1, 1.5, 2, 2.5, 2.5, 2, 1.5, 1]
            pod_percentile_labels = [f'PoD {p}th percentile' for p in percentiles]
            self._cache['pod_percentiles'] = (pod_percentiles, percentiles, pod_widths, pod_percentile_labels)
        return self._cache['pod_percentiles']

    def get_response_data(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get response data for plotting."""
        if 'response_data' not in self._cache:
            conc = np.array(self.df['conc'])
            conc_index = np.array(self.df['conc_index'])
            count = np.array(self.df[self.probe]['count'])
            total_count = np.array(self.df['total_count'])
            median_total_count = np.median(total_count)

            treatment_mask = conc_index > 0
            control_mask = conc_index == 0

            treatment_x = conc[conc_index[treatment_mask] - 1]
            treatment_y = (count[treatment_mask] / total_count[treatment_mask]) * median_total_count
            control_y = (count[control_mask] / total_count[control_mask]) * median_total_count
            control_y = filter_similar_control_lines(control_y, tolerance=0.02, min_lines=2)
            response_x = np.array(self.df[self.probe]['x'])
            response = np.array(self.df[self.probe]['response'])

            self._cache['response_data'] = (treatment_x, treatment_y, control_y, response_x, response)
        return self._cache['response_data']

def create_summary_table_data(probes: List[str], df: pd.Series, stats: Dict[str, np.ndarray],
                            weights: Dict[str, np.ndarray], conc_units: str, sort_by_abs_fc: bool = False) -> Dict[str, Dict[str, str]]:
    """Create summary table data for a list of probes."""
    # First create unsorted data
    unsorted_data = {}
    for probe in probes:
        mean_pod = np.mean(df[probe]['pod'])
        weight = weights['weight'][weights['probe'] == probe][0] if probe in weights['probe'] else 0.0
        l2fc = stats['l2fc'][stats['probe'] == probe][0]
        abs_fc = abs(l2fc)  # Calculate absolute fold change
        unsorted_data[probe] = {
            '_abs_fc': abs_fc,  # Keep as float for sorting
            'CDS': f"{df[probe]['cds']:.3f}",
            'Mean PoD': f"{10**mean_pod:.2g} {conc_units}",
            'Log2 Fold Change': f"{l2fc:.2f}",
            'Global PoD Weight': f"{weight:.3f}",
            'Response Range': f"{df[probe]['response_threshold_lower']:.1f} - {df[probe]['response_threshold_upper']:.1f}"
        }

    if sort_by_abs_fc:
        # Sort by absolute fold change in descending order
        sorted_probes = sorted(unsorted_data.keys(),
                             key=lambda x: unsorted_data[x]['_abs_fc'],
                             reverse=True)
        # Create new dictionary with sorted order
        table_data = {}
        for probe in sorted_probes:
            data = unsorted_data[probe].copy()
            data['_abs_fc'] = f"{data['_abs_fc']:.3f}"  # Convert to string after sorting
            table_data[probe] = data
    else:
        # If not sorting by abs_fc, just convert _abs_fc to string
        table_data = {probe: {k: (f"{v:.3f}" if k == '_abs_fc' else v)
                            for k, v in data.items()}
                     for probe, data in unsorted_data.items()}

    return table_data

def create_table_plot(data: Dict[str, Dict[str, str]], headers: Dict[str, Dict[str, Any]],
                     table_id: str, title: str, sort_by_abs_fc: bool = False) -> table.plot:
    """Create a MultiQC table plot with common configuration."""
    pconfig = {
        'id': table_id,
        'title': title,
        'namespace': 'BIFROST',
        'no_violin': True,
        'scale': False,
        'sort_rows': False,  # Disable automatic sorting
        'col1_header': 'Probe'
    }

    if sort_by_abs_fc:
        # Add _abs_fc to headers with supported options only
        headers['_abs_fc'] = {
            'title': '_abs_fc',
            'hidden': True,
            'description': 'Absolute fold change (for sorting)',
            'placement': 0  # Ensure it's the first column for sorting
        }

    return table.plot(data=data, headers=headers, pconfig=pconfig)

def create_probe_plot(df: pd.Series, probe: str, conc_units: str) -> linegraph.plot:
    """Creates a concentration-response plot for a specific probe."""
    logger.info(f"Creating concentration-response plot for probe {probe}")
    start_time = time.time()

    # Use ProbeData helper class
    probe_data = ProbeData(df, probe, conc_units)
    treatment_x, treatment_y, control_y, response_x, response = probe_data.get_response_data()

    # Calculate ymax
    ymax = float(max(np.max(treatment_y), np.max(control_y), np.max(response[2])) * 1.1)

    # Create plot data
    plot_data = {}
    extra_series = []

    # Add control lines
    extra_series.extend([
        {
            'name': 'Solvent control' if i == 0 else None,
            'pairs': [(float(10**(df['conc'][0] - 1)), float(y)), (float(10**(df['conc'][-1] + 1)), float(y))],
            'color': '#CCCCCC',
            'width': 1,
            'dash': 'dash',
            'showlegend': True if i == 0 else False
        }
        for i, y in enumerate(control_y)
    ])

    # Add PoD distribution lines if available
    pod_data = probe_data.pod_percentiles
    if pod_data is not None:
        pod_percentiles, _, pod_widths, pod_percentile_labels = pod_data
        extra_series.extend([
            {
                'name': pod_percentile_labels[i],
                'pairs': [(float(10**p), 0), (float(10**p), ymax)],
                'color': '#B19CD9',
                'width': pod_widths[i],
                'dash': 'solid',
                'showlegend': False
            }
            for i, p in enumerate(pod_percentiles)
        ])

    # Add response data
    extra_series.extend([
        {
            'name': '90% credible interval',
            'pairs': [(float(10**x), float(y)) for x, y in zip(response_x, response[0])],
            'color': '#FF8080',
            'width': 2,
            'dash': 'dash',
            'showlegend': True
        },
        {
            'name': None,
            'pairs': [(float(10**x), float(y)) for x, y in zip(response_x, response[2])],
            'color': '#FF8080',
            'width': 2,
            'dash': 'dash',
            'showlegend': False
        },
        {
            'name': 'Median response',
            'pairs': [(float(10**x), float(y)) for x, y in zip(response_x, response[1])],
            'color': '#FF0000',
            'width': 2,
            'showlegend': True
        },
        {
            'name': 'Treatment data',
            'pairs': [(float(10**x), float(y)) for x, y in zip(treatment_x, treatment_y)],
            'color': '#000000',
            'width': 0,
            'marker': 'x',
            'showlegend': True
        }
    ])

    # Add PoD mean line if available
    mean_pod = probe_data.mean_pod
    if mean_pod is not None:
        extra_series.append({
            'name': 'Mean PoD | Response',
            'pairs': [(float(10**mean_pod), 0), (float(10**mean_pod), ymax)],
            'color': '#663399',
            'width': 2,
            'dash': 'solid',
            'showlegend': True
        })

    # Create plot
    plot = linegraph.plot(
        plot_data,
        pconfig={
            'id': f'conc_response_{probe}',
            'title': f'{probe}',
            'xlab': f'Concentration ({conc_units})',
            'ylab': 'Normalised count',
            'xlog': True,
            'xmin': float(10**(df['conc'][0] - 1)),
            'ymin': 0,
            'ymax': ymax,
            'style': 'lines',
            'height': 400,
            'showlegend': True,
            'x_decimals': 2,
            'y_decimals': 0,
            'extra_series': extra_series
        }
    )

    elapsed_time = time.time() - start_time
    logger.info(f"Completed concentration-response plot for {probe} in {elapsed_time:.2f} seconds")
    return plot

def create_diagnostic_table_data(df: pd.Series, conc_units: str) -> Dict[str, Dict[str, Any]]:
    """Create diagnostic table data for all probes."""
    diagnostic_data = {}
    for probe in df['probes']:
        diag_text = df[probe]['diagnostics']
        probe_data = ProbeData(df, probe, conc_units)

        # Parse individual checks
        checks = {
            'Treedepth': '✓' if 'Treedepth satisfactory' in diag_text else '✗',
            'Divergences': '✓' if 'No divergent transitions' in diag_text else '✗',
            'E-BFMI': '✓' if 'E-BFMI satisfactory' in diag_text else '✗',
            'ESS': '✓' if 'effective sample size satisfactory' in diag_text else '✗',
            'R-hat': '✓' if 'R-hat greater than 1.01' not in diag_text else '✗'
        }

        # Calculate biological relevance score
        cds = probe_data.cds
        mean_pod = probe_data.mean_pod if probe_data.mean_pod is not None else float('inf')

        bio_score = 1 if cds > 0.5 else 0
        if not np.isinf(mean_pod):
            bio_score += (df['max_conc'] - mean_pod) / df['max_conc']

        # Extract R-hat parameters if present
        rhat_params = []
        if 'R-hat greater than 1.01' in diag_text:
            start_idx = diag_text.find('greater than 1.01:') + len('greater than 1.01:')
            end_idx = diag_text.find('Such high values')
            if start_idx > 0 and end_idx > start_idx:
                params_text = diag_text[start_idx:end_idx].strip()
                rhat_params = [p.strip() for p in params_text.split() if p.strip()]

        # Check for regularization recommendation
        needs_regularization = 'You should consider regularizating your model with additional prior information or a more effective parameterization' in diag_text

        # Format Mean PoD
        if not np.isnan(mean_pod):
            mean_pod_value = 10**mean_pod
            mean_pod_str = f"{mean_pod_value:.2g} {conc_units}"
        else:
            mean_pod_value = float('inf')
            mean_pod_str = "No response"

        # Add to diagnostic data
        diagnostic_data[probe] = {
            'CDS': float(cds),
            'CDS_str': f"{cds:.3f}",
            'Mean PoD': mean_pod_value,
            'Mean PoD_str': mean_pod_str,
            'Treedepth': checks['Treedepth'],
            'Divergences': checks['Divergences'],
            'E-BFMI': checks['E-BFMI'],
            'ESS': checks['ESS'],
            'R-hat': checks['R-hat'],
            'High R-hat Parameters': str(len(rhat_params)) if rhat_params else '0',
            'Response Range': f"{df[probe]['response_threshold_lower']:.1f} - {df[probe]['response_threshold_upper']:.1f}",
            'Needs Regularization': '⚠️' if needs_regularization else '✓',
            '_sort_score': bio_score
        }

    return diagnostic_data

@contextmanager
def timeout(seconds):
    """Context manager to enforce a timeout on a block of code."""
    def handler(signum, frame):
        raise TimeoutError(f"Operation timed out after {seconds} seconds")

    # Set the signal handler and a timer
    original_handler = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)

    try:
        yield
    finally:
        # Restore the original handler and cancel the alarm
        signal.alarm(0)
        signal.signal(signal.SIGALRM, original_handler)

def create_multiqc_report(summary_file, test_substance, cell_type, timepoint, conc_units, output_name,
                         interactive_plots=False, n_fold_change_probes=5, cds_threshold=0.5,
                         n_lowest_means=10, n_pod_stats=100, control_line_tolerance=0.02,
                         min_control_lines=2, plot_height=400, pod_vs_fc_height=600,
                         report_timeout=300, plots_force_flat_numseries=10000):
    """Create a MultiQC report from BIFROST data."""
    logger.info(f"Starting report generation for {test_substance} on {cell_type}")

    # Configure MultiQC for interactive plots if requested
    if interactive_plots:
        os.environ['MULTIQC_PLOTS_FORCE_INTERACTIVE'] = 'true'
        os.environ['MULTIQC_PLOTS_FLAT_NUMSERIES'] = str(plots_force_flat_numseries)

    # Load and process data
    logger.info(f"Loading summary data from {summary_file}")
    df = pd.read_json(
        summary_file,
        typ='series',
        orient='index',
        compression='zip',
        dtype_backend='numpy_nullable'
    )

    # Compute summary statistics and global PoD
    stats = calculate_summary_statistics(df)
    stats_filtered = filter_summary_statistics(stats, cds_threshold=cds_threshold)
    global_pod = get_global_pod(stats)
    weights = get_min_probe_weights(global_pod)

    # Initialize MultiQC
    multiqc.reset()
    multiqc.config.plots_force_flat = not interactive_plots
    multiqc.config.plots_flat_numseries = plots_force_flat_numseries if interactive_plots else 1000
    multiqc.config.skip_generalstats = True
    multiqc.config.skip_plots = False
    multiqc.config.skip_cleanup = True

    # Verify interactive plots configuration after initialization
    if interactive_plots:
        multiqc.config.plots_force_interactive = True
        multiqc.config.plots_flat_numseries = plots_force_flat_numseries

    # Create summary table
    logger.info("Creating summary table...")
    start_time = time.time()
    # Create summary table - format data as a dictionary of samples (metrics)
    summary_data = {
        'Global PoD': {
            'Value': f"{global_pod['global_pod']:.2g} {conc_units}"
        },
        'Maximum tested concentration': {
            'Value': f"{10**stats['max_conc']:.2g} {conc_units}"
        },
        'Number of probes analyzed': {
            'Value': str(len(stats['probe']))
        },
        'Number of hits': {
            'Value': str(int(global_pod['num_hits']))
        },
        'Number of CDS>0.5': {
            'Value': str(int(np.sum(stats['cds'] > 0.5)))
        },
        'Number of CDS=1.0': {
            'Value': str(int(np.sum(stats['cds'] == 1.0)))
        }
    }

    # Add minimum responding probe(s)
    valid_probes = weights['probe'][weights['probe'] != 'Max. conc.']
    if len(valid_probes) > 0:
        # Sort by weight and get the top probe
        order = np.argsort(weights['weight'][weights['probe'] != 'Max. conc.'])[::-1]
        top_probe = valid_probes[order][0]
        top_weight = weights['weight'][weights['probe'] != 'Max. conc.'][order][0]
        top_cds = weights['cds'][weights['probe'] != 'Max. conc.'][order][0]
        top_pod = weights['min_mean'][weights['probe'] != 'Max. conc.'][order][0]

        summary_data['Minimum responding probe'] = {
            'Value': f"{top_probe}, weight={top_weight:.2g}, CDS={top_cds:.2g}, Mean PoD={top_pod:.2g} {conc_units}"
        }

    # Add largest fold changes
    if len(stats['l2fc']) > 0:
        max_fc_idx = np.argmax(stats['l2fc'])
        min_fc_idx = np.argmin(stats['l2fc'])
        summary_data['Largest fold increase'] = {
            'Value': stats['probe'][max_fc_idx]
        }
        summary_data['Largest fold decrease'] = {
            'Value': stats['probe'][min_fc_idx]
        }

    # Add summary table to report
    summary_table = table.plot(
        data=summary_data,
        headers={
            'Value': {'title': 'Value'}
        },
        pconfig={
            'id': 'bifrost_summary',
            'title': 'BIFROST Analysis Summary',
            'namespace': 'BIFROST',
            'no_violin': True,
            'scale': False,  # Disable automatic scaling and coloring
            'sort_rows': False,
            'col1_header': 'Metric'  # This will label the first column as "Metric"
        }
    )

    # Create PoD vs fold-change plot data
    pod_vs_fc_data = {
        'Probe': [{  # Empty string as dataset name to avoid showing in tooltips
            'x': float(x),  # Convert numpy float to Python float
            'y': float(y),  # Convert numpy float to Python float
            'text': str(text),  # Convert numpy string to Python string
            'name': str(text)  # Add name for hover text
        } for x, y, text in zip(10**stats['pod'], stats['l2fc'], stats['probe'])]
    }

    # Calculate y-axis limits
    if len(stats['l2fc']) > 0:
        ymin, ymax = min(stats['l2fc'].min(), -2) - 1, max(stats['l2fc'].max(), 2) + 1
    else:
        ymin, ymax = -2, 2

    # Calculate x-axis maximum
    xmax = 10**np.array(stats['conc']).max() * 2  # Double the max concentration for padding

    # Create scatter plot using MultiQC's scatter plot type
    pod_vs_fc_plot = scatter.plot(
        pod_vs_fc_data,
        pconfig={
            'id': 'pod_vs_fc',
            'title': 'PoD vs Fold Change',
            'xlab': f'Mean PoD | Response ({conc_units})',
            'ylab': 'Max./min. log2 fold-change',
            'xlog': True,
            'xmin': 0,          # Set axis minimum to 0
            'xmax': xmax,       # Set axis maximum based on data
            'x_clipmin': 0.01,  # Clip data points below 0.01
            'x_clipmax': xmax,   # Clip data points above 100
            'x_decimals': 2,    # Format x-axis labels with 2 decimal places
            'ymin': ymin,       # Set y-axis minimum
            'ymax': ymax,       # Set y-axis maximum
            'marker_size': 5,
            'marker_line_width': 1,
            'color': 'black',   # Use color instead of marker_line_color
            'opacity': 1.0,     # Set full opacity
            'showlegend': False,  # Hide legend
            'height': pod_vs_fc_height,  # Make plot taller to accommodate labels
            'x_lines': [  # Add vertical lines using x_lines
                {
                    'value': float(global_pod['global_pod']),
                    'color': '#FF0000',  # Red
                    'width': 1,
                    'dash': 'solid',
                    'label': 'Global PoD'
                }
            ] + [
                {
                    'value': float(conc),
                    'color': '#D3D3D3',  # Light gray
                    'width': 1,
                    'dash': 'dash'
                }
                for conc in 10**np.array(stats['conc'])
            ],
            'y_lines': [  # Add horizontal line at y=0
                {
                    'value': 0,
                    'color': '#CCCCCC',  # Light gray
                    'width': 1,
                    'dash': 'dash'
                }
            ]
        }
    )

    # Create main BIFROST module for introduction and summary
    main_module = multiqc.BaseMultiqcModule(
        name='General',
        anchor='bifrost',
        href='https://github.com/your-repo/bifrost',
        info='BIFROST HTTr Analysis Report'
    )

    # Add introduction to main module
    main_module.add_section(
        name='Introduction',
        anchor='bifrost_intro',
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
        """
    )

    # Add summary section to main module
    main_module.add_section(
        name='Summary Statistics',
        anchor='bifrost_summary',
        plot=summary_table,
        description='''
        <p>Summary statistics from BIFROST analysis.</p>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li><strong>Global PoD</strong>: The estimated minimum effect concentration for any gene, summarizing the overall sensitivity of the system.</li>
            <li><strong>Maximum tested concentration</strong>: The highest concentration tested in the experiment.</li>
            <li><strong>Number of probes analyzed</strong>: Total number of probes included in the analysis.</li>
            <li><strong>Number of hits</strong>: Expected number of probes with a significant response.</li>
            <li><strong>Number of CDS>0.5 / CDS=1.0</strong>: Probes with strong concentration-dependent responses (CDS > 0.5) or maximal response (CDS = 1.0).</li>
            <li><strong>Minimum responding probe</strong>: The probe with the highest weight in the global PoD calculation, indicating the most sensitive response.</li>
            <li><strong>Largest fold increase/decrease</strong>: Probes with the largest positive or negative changes in expression.</li>
        </ul>
        '''
    )

    # Add PoD vs fold-change section to main module
    main_module.add_section(
        name='PoD vs Fold Change',
        anchor='bifrost_pod_vs_fc',
        plot=pod_vs_fc_plot,
        description='''
        <p>Maximum fold-change in expression over the tested concentration-range plotted against the
        probe-level PoD (mean given response). The red vertical line indicates the global PoD.
        Vertical grey lines are placed at the experimental test substance concentrations.</p>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li>Each point represents a probe.</li>
            <li><strong>X-axis (Mean PoD | Response)</strong>: Lower values indicate probes that respond at lower concentrations (more sensitive).</li>
            <li><strong>Y-axis (Max./min. log2 fold-change)</strong>: Higher absolute values indicate larger changes in expression.</li>
            <li>Probes in the lower-left region are most sensitive and show strong responses at low concentrations.</li>
            <li>The red vertical line (Global PoD) helps identify probes responding below the overall effect threshold.</li>
        </ul>
        '''
    )

    # Create module for probes with non-zero global PoD weight
    weighted_module = multiqc.BaseMultiqcModule(
        name='Probes with Non-zero Global PoD Weight',
        anchor='bifrost_weighted',
        info='Concentration-response plots for probes contributing to global PoD'
    )

    # Add description section to weighted module
    weighted_module.add_section(
        name='Overview',
        anchor='bifrost_weighted_overview',
        description='''
        <p>Concentration-response plots for probes that contribute to the global PoD calculation.</p>

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>These probes were selected based on their contribution to the global PoD calculation.</li>
            <li>Each probe's weight in the global PoD calculation is determined by its CDS (Concentration-Dependency Score) and its position in the PoD distribution.</li>
            <li>Probes are sorted by their weight in descending order, showing the most influential probes first.</li>
            <li>Only probes with non-zero weights are included, as these are the ones that meaningfully contribute to the global PoD estimate.</li>
        </ul>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li>These probes have the highest influence on the global PoD estimate.</li>
            <li>Their response curves and PoD distributions are most relevant for understanding the overall system sensitivity.</li>
            <li>The weight of each probe indicates its relative importance in determining the global PoD.</li>
        </ul>
        '''
    )

    # Add plots for probes with non-zero global PoD weight to weighted module
    valid_probes = weights['probe'][weights['probe'] != 'Max. conc.']
    probes_to_plot = valid_probes[np.argsort(weights['weight'][weights['probe'] != 'Max. conc.'])]
    logger.info(f"Found {len(probes_to_plot)} probes with non-zero global PoD weight to plot")

    if len(probes_to_plot) > 0:
        # Create summary table for weighted probes
        weighted_table_data = create_summary_table_data(probes_to_plot, df, stats, weights, conc_units, sort_by_abs_fc=True)
        weighted_summary_table = create_table_plot(
            data=weighted_table_data,
            headers={
                'CDS': {'title': 'CDS', 'description': 'Concentration-Dependency Score'},
                'Mean PoD': {'title': f'Mean PoD ({conc_units})', 'description': 'Mean point of departure'},
                'Log2 Fold Change': {'title': 'Log2 Fold Change', 'description': 'Maximum fold change in expression'},
                'Global PoD Weight': {'title': 'Global PoD Weight', 'description': 'Weight in global PoD calculation'},
                'Response Range': {'title': 'Response Range', 'description': 'Range of response thresholds'}
            },
            table_id='bifrost_weighted_summary',
            title='Summary Statistics for Probes with Non-zero Global PoD Weight',
            sort_by_abs_fc=True
        )

        # Add summary table to weighted module
        weighted_module.add_section(
            name='Probe Summary Statistics',
            anchor='bifrost_weighted_summary',
            plot=weighted_summary_table,
            description='''
            <p>Summary statistics for probes contributing to the global PoD calculation.</p>
            <p>The table shows key metrics for each probe, sorted by their weight in the global PoD calculation.</p>
            '''
        )

        logger.info("Generating concentration-response plots for probes with non-zero global PoD weight...")
        start_time = time.time()
        for i, probe in enumerate(probes_to_plot, 1):
            logger.info(f"Plotting probe {i}/{len(probes_to_plot)}: {probe}")
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            weighted_module.add_section(
                name=probe,
                anchor=f'bifrost_weighted_{probe}',
                plot=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}'
            )
        logger.info(f"Completed all concentration-response plots for non-zero global PoD weight probes in {time.time() - start_time:.2f} seconds")

    # Create module for probes with largest fold changes
    fc_module = multiqc.BaseMultiqcModule(
        name='Probes with Largest Fold Changes',
        anchor='bifrost_fc',
        info='Concentration-response plots for probes with extreme expression changes'
    )

    # Add description section to fc module
    fc_module.add_section(
        name='Overview',
        anchor='bifrost_fc_overview',
        description=f'''
        <p>Concentration-response plots for probes with the most extreme expression changes, separated into upregulated and downregulated genes.</p>

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>This module shows the probes with the most extreme fold changes in expression, divided into two categories:
                <ul>
                    <li><strong>Most Upregulated</strong>: The {n_fold_change_probes} probes with the largest positive fold changes (increased expression)</li>
                    <li><strong>Most Downregulated</strong>: The {n_fold_change_probes} probes with the largest negative fold changes (decreased expression)</li>
                </ul>
            </li>
            <li>Fold changes are calculated as log2 ratios of expression at the maximum response concentration compared to the control.</li>
            <li>Probes are selected regardless of their CDS or PoD values, focusing solely on the magnitude of expression change.</li>
        </ul>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li>These probes show the most extreme changes in expression, regardless of sensitivity.</li>
            <li>Useful for identifying outliers or highly dynamic responses.</li>
            <li>Note that large fold changes don't necessarily indicate biological relevance - check the CDS and PoD values for context.</li>
        </ul>
        '''
    )

    # Add section for most upregulated probes
    fc_module.add_section(
        name='Most Upregulated Probes',
        anchor='bifrost_fc_up',
        description=f'''
        <p>Concentration-response plots for the {n_fold_change_probes} probes with the largest positive fold changes (increased expression).</p>

        <p><strong>Selection Details:</strong></p>
        <ul>
            <li>These probes show the strongest increase in expression across the concentration range.</li>
            <li>Selected based on the maximum positive log2 fold change relative to control.</li>
            <li>Sorted by fold change magnitude in descending order.</li>
        </ul>
        '''
    )

    # Add section for most downregulated probes
    fc_module.add_section(
        name='Most Downregulated Probes',
        anchor='bifrost_fc_down',
        description=f'''
        <p>Concentration-response plots for the {n_fold_change_probes} probes with the largest negative fold changes (decreased expression).</p>

        <p><strong>Selection Details:</strong></p>
        <ul>
            <li>These probes show the strongest decrease in expression across the concentration range.</li>
            <li>Selected based on the maximum negative log2 fold change relative to control.</li>
            <li>Sorted by fold change magnitude in ascending order (most negative first).</li>
        </ul>
        '''
    )

    # Add plots for most upregulated probes
    if len(stats['l2fc']) > 0:
        # Sort by absolute fold change magnitude
        abs_fc = np.abs(stats['l2fc'])
        index = np.argsort(abs_fc)[::-1]  # Sort in descending order
        n_up = min(n_fold_change_probes, len(stats['l2fc']))
        # Get probes with largest absolute fold changes that are positive
        up_probes = stats['probe'][index][stats['l2fc'][index] > 0][:n_up]
        logger.info(f"Found {len(up_probes)} probes with largest positive fold changes to plot")

        # Create summary table for upregulated probes
        up_table_data = create_summary_table_data(up_probes, df, stats, weights, conc_units, sort_by_abs_fc=True)
        up_summary_table = create_table_plot(
            data=up_table_data,
            headers={
                'CDS': {'title': 'CDS', 'description': 'Concentration-Dependency Score'},
                'Mean PoD': {'title': f'Mean PoD ({conc_units})', 'description': 'Mean point of departure'},
                'Log2 Fold Change': {'title': 'Log2 Fold Change', 'description': 'Maximum positive fold change in expression'},
                'Global PoD Weight': {'title': 'Global PoD Weight', 'description': 'Weight in global PoD calculation'},
                'Response Range': {'title': 'Response Range', 'description': 'Range of response thresholds'}
            },
            table_id='bifrost_fc_up_summary',
            title='Summary Statistics for Most Upregulated Probes',
            sort_by_abs_fc=True
        )

        # Add summary table to upregulated section
        fc_module.add_section(
            name='Upregulated Probes Summary',
            anchor='bifrost_fc_up_summary',
            plot=up_summary_table,
            description='''
            <p>Summary statistics for the most upregulated probes.</p>
            <p>The table shows key metrics for each probe, sorted by their fold change magnitude.</p>
            '''
        )

        logger.info("Generating concentration-response plots for most upregulated probes...")
        start_time = time.time()
        for i, probe in enumerate(up_probes, 1):
            logger.info(f"Plotting probe {i}/{len(up_probes)}: {probe}")
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            fc_module.add_section(
                name=probe,
                anchor=f'bifrost_fc_up_{probe}',
                plot=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}, Log2 Fold Change = {stats["l2fc"][stats["probe"] == probe][0]:.2f}'
            )
        logger.info(f"Completed all concentration-response plots for most upregulated probes in {time.time() - start_time:.2f} seconds")

    # Add plots for most downregulated probes
    if len(stats['l2fc']) > 0:
        # Sort by absolute fold change magnitude
        abs_fc = np.abs(stats['l2fc'])
        index = np.argsort(abs_fc)[::-1]  # Sort in descending order
        n_down = min(n_fold_change_probes, len(stats['l2fc']))
        # Get probes with largest absolute fold changes that are negative
        down_probes = stats['probe'][index][stats['l2fc'][index] < 0][:n_down]
        logger.info(f"Found {len(down_probes)} probes with largest negative fold changes to plot")

        # Create summary table for downregulated probes
        down_table_data = create_summary_table_data(down_probes, df, stats, weights, conc_units, sort_by_abs_fc=True)
        down_summary_table = create_table_plot(
            data=down_table_data,
            headers={
                'CDS': {'title': 'CDS', 'description': 'Concentration-Dependency Score'},
                'Mean PoD': {'title': f'Mean PoD ({conc_units})', 'description': 'Mean point of departure'},
                'Log2 Fold Change': {'title': 'Log2 Fold Change', 'description': 'Maximum negative fold change in expression'},
                'Global PoD Weight': {'title': 'Global PoD Weight', 'description': 'Weight in global PoD calculation'},
                'Response Range': {'title': 'Response Range', 'description': 'Range of response thresholds'}
            },
            table_id='bifrost_fc_down_summary',
            title='Summary Statistics for Most Downregulated Probes',
            sort_by_abs_fc=True
        )

        # Add summary table to downregulated section
        fc_module.add_section(
            name='Downregulated Probes Summary',
            anchor='bifrost_fc_down_summary',
            plot=down_summary_table,
            description='''
            <p>Summary statistics for the most downregulated probes.</p>
            <p>The table shows key metrics for each probe, sorted by their fold change magnitude.</p>
            '''
        )

        logger.info("Generating concentration-response plots for most downregulated probes...")
        start_time = time.time()
        for i, probe in enumerate(down_probes, 1):
            logger.info(f"Plotting probe {i}/{len(down_probes)}: {probe}")
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            fc_module.add_section(
                name=probe,
                anchor=f'bifrost_fc_down_{probe}',
                plot=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}, Log2 Fold Change = {stats["l2fc"][stats["probe"] == probe][0]:.2f}'
            )
        logger.info(f"Completed all concentration-response plots for most downregulated probes in {time.time() - start_time:.2f} seconds")

    # Create module for lowest means with CDS > 0.5
    lowest_means_module = multiqc.BaseMultiqcModule(
        name='Lowest Mean PoDs (CDS > 0.5)',
        anchor='bifrost_lowest_means',
        info='Concentration-response plots for most sensitive probes'
    )

    # Add description section to lowest means module
    lowest_means_module.add_section(
        name='Overview',
        anchor='bifrost_lowest_means_overview',
        description=f'''
        <p>Concentration-response plots for the most sensitive probes with strong evidence of response.</p>

        <p><strong>Probe Selection:</strong></p>
        <ul>
            <li>This section displays the {n_lowest_means} probes with the lowest mean PoDs (most sensitive) that meet two criteria:
                <ul>
                    <li>CDS > {cds_threshold} (strong evidence for a concentration-dependent response)</li>
                    <li>Valid PoD estimate (mean PoD less than maximum tested concentration)</li>
                </ul>
            </li>
            <li>Probes are first filtered to include only those with CDS > {cds_threshold}, ensuring reliable concentration-dependent responses.</li>
            <li>Among these filtered probes, the {n_lowest_means} with the lowest mean PoDs are selected.</li>
            <li>If fewer than {n_lowest_means} probes meet these criteria, all qualifying probes are shown.</li>
        </ul>

        <p><strong>How to interpret:</strong></p>
        <ul>
            <li>These are the most sensitive probes (lowest mean PoD) with strong evidence for a response (CDS > {cds_threshold}).</li>
            <li>Useful for identifying the earliest responding genes.</li>
            <li>The combination of low PoD and high CDS suggests these are reliable early indicators of biological response.</li>
        </ul>
        '''
    )

    # Add plots for lowest means to lowest means module
    n_probe = len(stats['probe'])
    probes_to_plot = stats['probe'][np.argsort(stats['pod'])][:min(n_probe, n_lowest_means)]
    logger.info(f"Found {len(probes_to_plot)} probes with lowest means and CDS > 0.5 to plot")

    if len(probes_to_plot) > 0:
        # Create summary table for lowest means probes
        lowest_means_table_data = create_summary_table_data(probes_to_plot, df, stats, weights, conc_units, sort_by_abs_fc=True)
        lowest_means_summary_table = create_table_plot(
            data=lowest_means_table_data,
            headers={
                'CDS': {'title': 'CDS', 'description': 'Concentration-Dependency Score'},
                'Mean PoD': {'title': f'Mean PoD ({conc_units})', 'description': 'Mean point of departure'},
                'Log2 Fold Change': {'title': 'Log2 Fold Change', 'description': 'Maximum fold change in expression'},
                'Global PoD Weight': {'title': 'Global PoD Weight', 'description': 'Weight in global PoD calculation'},
                'Response Range': {'title': 'Response Range', 'description': 'Range of response thresholds'}
            },
            table_id='bifrost_lowest_means_summary',
            title='Summary Statistics for Most Sensitive Probes (CDS > 0.5)',
            sort_by_abs_fc=True
        )

        # Add summary table to lowest means module
        lowest_means_module.add_section(
            name='Probe Summary Statistics',
            anchor='bifrost_lowest_means_summary',
            plot=lowest_means_summary_table,
            description='''
            <p>Summary statistics for the most sensitive probes with CDS > 0.5.</p>
            <p>The table shows key metrics for each probe, sorted by their mean PoD (most sensitive first).</p>
            '''
        )

        logger.info("Generating concentration-response plots for probes with lowest means...")
        start_time = time.time()
        for i, probe in enumerate(probes_to_plot, 1):
            logger.info(f"Plotting probe {i}/{len(probes_to_plot)}: {probe}")
            conc_response_plot = create_probe_plot(df, probe, conc_units)
            lowest_means_module.add_section(
                name=probe,
                anchor=f'bifrost_lowest_means_{probe}',
                plot=conc_response_plot,
                description=f'CDS = {df[probe]["cds"]:.3f}, Mean PoD = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}'
            )
        logger.info(f"Completed all concentration-response plots for lowest means probes in {time.time() - start_time:.2f} seconds")

    # Create module for PoD statistics
    stats_module = multiqc.BaseMultiqcModule(
        name='Probe-level PoD Statistics',
        anchor='bifrost_stats',
        info=f'Detailed statistics for probes with CDS > {cds_threshold}'
    )

    # Add PoD statistics table to stats module
    if n_probe > 0:
        # Get probes with CDS > 0.5 and sort by PoD
        cds_mask = stats['cds'] > cds_threshold
        probes = stats['probe'][cds_mask][np.argsort(stats['pod'][cds_mask])][:n_pod_stats]
        logger.info(f"Found {len(probes)} probes with CDS > 0.5 to include in PoD statistics table (limited to top {n_pod_stats})")
        cds = stats['cds'][cds_mask][np.argsort(stats['pod'][cds_mask])][:n_pod_stats]

        # Create table data
        table_data = {}
        for i, (probe, cds_val) in enumerate(zip(probes, cds)):
            n_pod_samples = len(df[probe]['pod'])
            extended_pod_samples = np.concatenate((df[probe]['pod'],
                                [df['max_conc'] for _ in range(df['n_samp'] - n_pod_samples)]))
            pod_percentiles = np.percentile(extended_pod_samples, q=(5, 25, 50, 75, 95))

            # Format PoD values
            pod_values = []
            for pod_val in pod_percentiles:
                if pod_val < df['max_conc']:
                    # Convert to integer and format without scientific notation
                    pod_values.append(f"{int(10**pod_val)}")
                else:
                    pod_values.append(f">{int(10**pod_val)}")

            # Add row to table data
            table_data[probe] = {
                'CDS': f"{cds_val:.3f}",
                '5th percentile': pod_values[0],
                '25th percentile': pod_values[1],
                '50th percentile': pod_values[2],
                '75th percentile': pod_values[3],
                '95th percentile': pod_values[4]
            }

        # Create table plot
        pod_stats_table = table.plot(
            data=table_data,
            headers={
                'CDS': {'title': 'CDS', 'format': '{:.3f}', 'description': 'Concentration-Dependency Score'},
                '5th percentile': {'title': f'5th percentile ({conc_units})', 'description': '5th percentile of PoD distribution'},
                '25th percentile': {'title': f'25th percentile ({conc_units})', 'description': '25th percentile of PoD distribution'},
                '50th percentile': {'title': f'50th percentile ({conc_units})', 'description': 'Median of PoD distribution'},
                '75th percentile': {'title': f'75th percentile ({conc_units})', 'description': '75th percentile of PoD distribution'},
                '95th percentile': {'title': f'95th percentile ({conc_units})', 'description': '95th percentile of PoD distribution'}
            },
            pconfig={
                'id': 'bifrost_stats_table',
                'title': 'Probe-level PoD Statistics (CDS > 0.5)',
                'namespace': 'BIFROST',
                'no_violin': True,
                'scale': False,  # Disable automatic scaling and coloring
                'sort_rows': False,
                'col1_header': 'Probe'  # Label first column as Probe
            }
        )

        # Add table to section
        stats_module.add_section(
            name='PoD Statistics Table',
            anchor='bifrost_stats_table',
            plot=pod_stats_table,
            description="""
            <p>Summary statistics for probes with CDS > 0.5, showing the distribution of PoD (Point of Departure) values.</p>
            <p>The table includes:</p>
            <ul>
                <li><strong>CDS</strong>: Concentration-Dependency Score (probability of response below max concentration)</li>
                <li><strong>PoD percentiles</strong>: Different quantiles of the PoD distribution</li>
                <li>Values are shown in {conc_units}</li>
                <li>Probes are sorted by median PoD (50th percentile)</li>
                <li>Only shows the top {n_pod_stats} probes with CDS > 0.5</li>
            </ul>

            <p>Hover over column headers for more detailed descriptions.</p>
            """
        )

    # Create diagnostic table data with parsed checks
    diagnostic_data = create_diagnostic_table_data(df, conc_units)

    # Create diagnostic table
    diagnostic_table = create_table_plot(
        data={k: {sk: v[sk] for sk in ['CDS_str', 'Mean PoD_str', 'Treedepth', 'Divergences', 'E-BFMI', 'ESS', 'R-hat', 'High R-hat Parameters', 'Response Range', 'Needs Regularization', '_sort_score']} for k, v in diagnostic_data.items()},
        headers={
            'CDS_str': {
                'title': 'CDS',
                'format': '{:.3f}',
                'description': f'Concentration-Dependency Score (probability of response below max concentration, threshold = {cds_threshold})',
                'cond_formatting_rules': {
                    'pass': [{'gt': cds_threshold}]  # Highlight probes with CDS > threshold
                }
            },
            'Mean PoD_str': {
                'title': f'Mean PoD ({conc_units})',
                'description': 'Mean point of departure (effect concentration). "No response" indicates no valid PoD samples.',
                'cond_formatting_rules': {
                    'warn': [{'s_eq': 'No response'}]  # Highlight probes with no response
                }
            },
            'Treedepth': {
                'title': 'Treedepth',
                'description': 'Sampler transitions treedepth check',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for pass
                    'fail': [{'s_eq': '✗'}]   # Red for fail
                }
            },
            'Divergences': {
                'title': 'Divergences',
                'description': 'Check for divergent transitions',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for pass
                    'fail': [{'s_eq': '✗'}]   # Red for fail
                }
            },
            'E-BFMI': {
                'title': 'E-BFMI',
                'description': 'HMC potential energy check',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for pass
                    'fail': [{'s_eq': '✗'}]   # Red for fail
                }
            },
            'ESS': {
                'title': 'ESS',
                'description': 'Effective sample size check',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for pass
                    'fail': [{'s_eq': '✗'}]   # Red for fail
                }
            },
            'R-hat': {
                'title': 'R-hat',
                'description': 'Gelman-Rubin convergence diagnostic',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for pass
                    'fail': [{'s_eq': '✗'}]   # Red for fail
                }
            },
            'High R-hat Parameters': {'title': '# Parameters with R-hat > 1.01', 'description': 'Number of parameters with high R-hat values'},
            'Response Range': {'title': 'Response Range', 'description': 'Range of response thresholds'},
            'Needs Regularization': {
                'title': '⚠️ Regularization',
                'description': 'Model may need regularization',
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '✓'}],  # Green for no regularization needed
                    'warn': [{'s_eq': '⚠️'}]  # Orange for needs regularization
                }
            },
            '_sort_score': {
                'title': '_sort_score',
                'hidden': True  # Hide the sorting column from display
            }
        },
        table_id='bifrost_diagnostics_table',
        title='Probe Diagnostic Summary'
    )

    # Create module for diagnostics
    diag_module = multiqc.BaseMultiqcModule(
        name='Diagnostic Summary',
        anchor='bifrost_diagnostics',
        info='Model diagnostics and quality checks'
    )

    # Add diagnostic table to diag module
    diag_module.add_section(
        name='Diagnostic Table',
        anchor='bifrost_diagnostics_table',
        plot=diagnostic_table,
        description=f'''
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
            <li><strong>R-hat</strong>: Checks Gelman-Rubin convergence diagnostic (should be < 1.01)</li>
        </ul>

        <p><strong>Additional Information:</strong></p>
        <ul>
            <li><strong>High R-hat Parameters</strong>: Lists parameters with R-hat > 1.01 that may need attention</li>
            <li><strong>CDS</strong>: Concentration-Dependency Score (probability of response below max concentration)</li>
            <li><strong>Mean PoD</strong>: Mean point of departure (effect concentration)</li>
            <li><strong>Response Range</strong>: Range of response thresholds</li>
            <li><strong>Regularization</strong>: ⚠️ indicates model may need regularization</li>
        </ul>

        <p>Probes with "No response" in the Mean PoD column did not show a significant response in the tested range.
        Failed diagnostic checks (red ✗) indicate potential model issues for that probe.</p>
        '''
    )

    # Add all modules to report
    multiqc.report.modules.extend([
        main_module,
        weighted_module,
        fc_module,
        lowest_means_module,
        stats_module,
        diag_module
    ])

    # Write report
    logger.info("Generating report (this may take a few minutes for large datasets)...")
    start_time = time.time()

    # Create a thread to monitor progress
    stop_monitoring = threading.Event()

    def monitor_progress():
        last_update = time.time()
        while not stop_monitoring.is_set():
            current_time = time.time()
            if current_time - last_update >= 30:  # Log every 30 seconds instead of 5
                elapsed = current_time - start_time
                logger.info(f"Still generating report... ({elapsed:.0f} seconds elapsed)")
                last_update = current_time
            time.sleep(1)

    # Start the monitoring thread
    monitor_thread = threading.Thread(target=monitor_progress, daemon=True)
    monitor_thread.start()

    try:
        multiqc.config.verbose = True
        try:
            with timeout(report_timeout):
                multiqc.write_report(
                    output_dir=os.path.dirname(output_name),
                    filename=os.path.basename(output_name),
                    title=f'BIFROST HTTr Analysis - {test_substance} ({cell_type})',
                    report_comment=f'Analysis of {test_substance} on {cell_type} cells after {timepoint} exposure',
                    force=True
                )
        except TimeoutError:
            logger.error(f"Report generation timed out after {report_timeout} seconds")
            logger.error("Try running with --interactive-plots option for faster rendering")
            sys.exit(1)

    except Exception as e:
        logger.error(f"Error during report generation: {str(e)}")
        raise
    finally:
        stop_monitoring.set()
        monitor_thread.join(timeout=1)

    elapsed_time = time.time() - start_time
    logger.info(f"Report generation complete in {elapsed_time:.0f} seconds")
    if elapsed_time > 60 and not interactive_plots:
        logger.info("Note: Consider using --interactive-plots for faster rendering with large datasets")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Create Bifrost HTTR reports using MultiQC')
    parser.add_argument('--summary-file',
                      required=True,
                      help='Path to the summary JSON file')
    parser.add_argument('--test-substance',
                      default='MyChemical',
                      help='Name of the test substance (default: MyChemical)')
    parser.add_argument('--cell-type',
                      default='MyCell',
                      help='Type of cell used in the test (default: MyCell)')
    parser.add_argument('--output-name',
                      default='multiqc_report.html',
                      help='Name for the output report (default: multiqc_report.html)')
    parser.add_argument('--timepoint',
                      default='24 hours',
                      help='Exposure duration within experiment (default: 24 hours)')
    parser.add_argument('--conc-units',
                      default='uM',
                      choices=['uM', 'ugml-1', 'mgml-1'],
                      help='Concentration units (default: uM)')
    parser.add_argument('--interactive-plots',
                      action='store_true',
                      help='Force interactive plots (may be faster for large datasets)')
    parser.add_argument('--n-fold-change-probes',
                      type=int,
                      default=5,
                      help='Number of most up/down regulated probes to show (default: 5)')

    # Add new parameters for configurable settings
    parser.add_argument('--cds-threshold',
                      type=float,
                      default=0.5,
                      help='Concentration-Dependency Score threshold for filtering probes (default: 0.5)')
    parser.add_argument('--n-lowest-means',
                      type=int,
                      default=10,
                      help='Number of lowest mean PoD probes to show (default: 10)')
    parser.add_argument('--n-pod-stats',
                      type=int,
                      default=100,
                      help='Number of probes to include in PoD statistics table (default: 100)')
    parser.add_argument('--control-line-tolerance',
                      type=float,
                      default=0.02,
                      help='Tolerance for filtering similar control lines (default: 0.02)')
    parser.add_argument('--min-control-lines',
                      type=int,
                      default=2,
                      help='Minimum number of control lines to show (default: 2)')
    parser.add_argument('--plot-height',
                      type=int,
                      default=400,
                      help='Height of concentration-response plots in pixels (default: 400)')
    parser.add_argument('--pod-vs-fc-height',
                      type=int,
                      default=600,
                      help='Height of PoD vs Fold Change plot in pixels (default: 600)')
    parser.add_argument('--report-timeout',
                      type=int,
                      default=300,
                      help='Timeout in seconds for report generation (default: 300)')
    parser.add_argument('--plots-force-flat-numseries',
                      type=int,
                      default=10000,
                      help='Maximum number of series for flat plots (default: 10000)')

    return parser.parse_args()

def main():
    """Main entry point for report generation script."""
    args = parse_args()

    try:
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
            args.report_timeout,
            args.plots_force_flat_numseries
        )
    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        raise

if __name__ == '__main__':
    main()
