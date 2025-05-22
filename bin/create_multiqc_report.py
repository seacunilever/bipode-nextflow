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
from multiqc.plots import table, linegraph, scatter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def calculate_summary_statistics(df: pd.Series) -> Dict[str, Union[np.ndarray, float, int]]:
    """Calculates summary statistics for BIFROST analysis.

    This function computes various summary statistics including PoD mean, log2 fold-change
    extrema, and CDS scores for each probe in the dataset.

    Args:
        df: BIFROST summary data as a pandas Series containing probe information.

    Returns:
        A dictionary containing:
            - probe: Array of probe identifiers
            - pod: Array of PoD means for each probe
            - cds: Array of CDS scores for each probe
            - l2fc: Array of log2 fold changes for each probe
            - max_conc: Maximum tested concentration
            - n_samp: Number of samples
            - conc: Array of concentration values

    Raises:
        KeyError: If required keys are missing from the input Series.
        ValueError: If probe data is malformed.
    """
    logger.info("Starting summary statistics calculation...")

    # Extend PoD sample arrays using maximum concentration value to redefine PoD percentile estimates
    probes = np.array(df['probes'])
    logger.info(f"Processing {len(probes)} probes...")

    pod = np.array([np.mean(df[i]['pod']) if len(df[i]['pod']) > 0 else np.nan for i in df['probes']])
    logger.info(f"Calculated PoD means for {np.sum(~np.isnan(pod))} probes")

    cds = np.array([df[i]['cds'] for i in probes])
    logger.info(f"Calculated CDS scores for {len(cds)} probes")

    l2fc = np.empty(probes.shape[0], dtype='float')
    logger.info("Calculating log2 fold changes...")
    for i, probe in enumerate(probes):
        y = np.array(df[probe]['response'][1])
        index = np.argmax(np.abs(np.log2(y / y[0])))
        l2fc[i] = np.log2(y[index] / y[0])
    logger.info("Log2 fold changes calculation complete")

    stats = {
        'probe': probes,
        'pod': pod,
        'cds': cds,
        'l2fc': l2fc,
        'max_conc': df['max_conc'],
        'n_samp': df['n_samp'],
        'conc': df['conc']
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

def create_multiqc_report(summary_file, test_substance, cell_type, timepoint, conc_units, output_name):
    """Create a MultiQC report from BIFROST data."""
    logger.info(f"Starting report generation for {test_substance} on {cell_type}")

    # Load and process data
    df = pd.read_json(summary_file, typ='series', orient='index', compression='zip')

    # Compute summary statistics and global PoD
    logger.info("Calculating summary statistics...")
    stats = calculate_summary_statistics(df)
    logger.info(f"Found {len(stats['probe'])} probes in summary statistics")

    logger.info("Filtering summary statistics...")
    stats_filtered = filter_summary_statistics(stats, cds_threshold=0.5)
    logger.info(f"After filtering: {len(stats_filtered['probe'])} probes")

    logger.info("Calculating global PoD...")
    global_pod = get_global_pod(stats)
    logger.info(f"Global PoD calculated: {global_pod['global_pod']}")

    logger.info("Calculating probe weights...")
    weights = get_min_probe_weights(global_pod)
    logger.info(f"Calculated weights for {len(weights['probe'])} probes")

    # Initialize MultiQC
    multiqc.reset()

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
            'height': 600,  # Make plot taller to accommodate labels
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

    # Add sections to report
    module = multiqc.BaseMultiqcModule(
        name='BIFROST',
        anchor='bifrost',
        href='https://github.com/your-repo/bifrost',
        info='BIFROST HTTr Analysis Report'
    )

    # Add introduction
    module.add_section(
        name='Introduction',
        anchor='bifrost_intro',
        description=f"""
        This report contains analysis of high-throughput transcriptomics data (HTTr) obtained after
        exposing {cell_type} cells for {timepoint} to {test_substance}. The BIFROST model
        (Bayesian inference for region of signal threshold) is a statistical model for analysis of
        HTTr concentration-response data.
        """
    )

    # Add summary section
    module.add_section(
        name='Summary Statistics',
        anchor='bifrost_summary',
        plot=summary_table,
        description='Summary statistics from BIFROST analysis.'
    )

    # Add PoD vs fold-change section
    module.add_section(
        name='PoD vs Fold Change',
        anchor='bifrost_pod_vs_fc',
        plot=pod_vs_fc_plot,
        description='''
        Maximum fold-change in expression over the tested concentration-range plotted against the
        probe-level PoD (mean given response). The red vertical line indicates the global PoD.
        '''
    )

    # Add module to report
    multiqc.report.modules.append(module)

    # Write report
    multiqc.write_report(
        output_dir=os.path.dirname(output_name),
        filename=os.path.basename(output_name),
        title=f'BIFROST HTTr Analysis - {test_substance} ({cell_type})',
        report_comment=f'Analysis of {test_substance} on {cell_type} cells after {timepoint} exposure',
        force=True  # Add force flag to overwrite existing reports
    )

    logger.info("Report generation complete")

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
            args.output_name
        )
    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        raise

if __name__ == '__main__':
    main()
