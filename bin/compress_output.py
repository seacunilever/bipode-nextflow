#!/usr/bin/env python

import argparse
import os
from pathlib import Path
from typing import Dict, Any, Union, List, Optional, Tuple
import numpy as np
import pandas as pd
import pickle
import subprocess


def get_confidence_threshold_probability_density(x: np.ndarray) -> np.ndarray:
    """
    Evaluate the probability density for the function describing uncertainty in CDS threshold.

    Args:
        x: Array of values at which to calculate density

    Returns:
        Array of corresponding probability densities
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


def get_minimum_pod_means(
    pod_means: np.ndarray,
    cds: np.ndarray,
    cds_thresholds: np.ndarray,
    max_conc: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute minimum PoD means for each value in the specified array.

    Args:
        pod_means: Array of PoD means
        cds: Array of CDS values
        cds_thresholds: Array of thresholds on which to filter probes
        max_conc: Maximum concentration value

    Returns:
        Tuple containing:
            - Array of minimum mean values
            - Array of corresponding probe IDs
            - Array of corresponding CDS values
    """
    min_means = np.full(len(cds_thresholds), max_conc)
    min_probes = np.full(len(cds_thresholds), 'Max. conc.', dtype='object')
    min_cds = np.full(len(cds_thresholds), 0, dtype='float')
    for i, j in enumerate(cds_thresholds):
        mask = cds >= j
        if np.sum(mask) > 0:
            pod = pod_means[mask]
            index = np.argmin(pod)
            min_means[i] = pod[index]
            min_probes[i] = cds[mask][index]
            min_cds[i] = cds[mask][index]

    return min_means, min_probes, min_cds


def get_global_pod(df: pd.Series, seed: Optional[int] = None) -> Dict[str, Any]:
    """
    Calculate and return global PoD.

    Args:
        df: BIFROST summary Series
        seed: Optional random seed for reproducibility

    Returns:
        Dictionary of global PoD-related statistics
    """
    if seed is not None:
        np.random.seed(seed)

    # Extract PoD means and CDS
    pod_means = np.array([np.mean(df[i]['pod']) if len(df[i]['pod']) > 0 else np.nan for i in df['probes']])
    cds = np.array([df[i]['cds'] for i in df['probes']])

    # Get dictionary of minimum probes
    dq = 0.025
    quantiles = np.arange(0.5, 1 + dq, dq)
    min_means, min_probes, min_cds = get_minimum_pod_means(pod_means, cds, quantiles, df['max_conc'])

    # Calculate weights and global PoD
    weights = get_confidence_threshold_probability_density(quantiles)
    weight_sum = np.sum(weights)
    global_pod = 10 ** (np.sum(min_means * weights) / weight_sum)

    # Calculate number of hits at each confidence threshold
    num_hits = np.array([np.sum(cds >= i) for i in quantiles])
    expected_num_hits = np.round(np.sum(num_hits * weights) / weight_sum)

    results = {'global_pod': global_pod, 'num_hits': expected_num_hits,
               'means': min_means, 'probes': min_probes, 'weights': weights,
               'quantiles': quantiles, 'cds': min_cds}

    return results


def compress_output(analysis_dir: Union[str, Path], path_to_summary: Union[str, Path], seed: Optional[int] = None, no_compression: bool = False) -> None:
    """
    Compress intermediate output into a single pandas DataFrame.

    Args:
        analysis_dir: Path to parent directory for intermediate output
        path_to_summary: Path to summary file
        seed: Optional random seed for reproducibility
        no_compression: If True, save as plain JSON without compression
    """
    if seed is not None:
        np.random.seed(seed)

    analysis_path = Path(analysis_dir)

    # Determine probe IDs
    data_files = [f for f in analysis_path.iterdir() if f.is_file()]
    probes = np.array([file.stem for file in data_files])

    # Create empty pandas series
    summary = pd.Series(dtype='object')

    # Extract details inputs universal to all chemicals/probes
    with open(analysis_path / f"{probes[0]}.pkl", 'rb') as f:
        data = pickle.load(f)

    for key in ['n_samp',
                'n_sample', 'n_treatment_batch', 'total_count', 'n_batch', 'batch_index',
                'n_conc', 'conc', 'conc_index', 'max_conc',
                ]:
        summary[key] = data[key]

    summary['probes'] = probes

    # Extract probe-specific information
    for probe in probes:
        with open(analysis_path / f"{probe}.pkl", 'rb') as f:
            data = pickle.load(f)

        summary[probe] = pd.Series(dtype='object')
        summary[probe]['diagnostics'] = data['diagnostics']
        summary[probe]['parameters'] = data['parameters']

        for par in data['fit'].index:
            summary[probe][par] = data['fit'][par]
        summary[probe]['count'] = data['count']

    # Calculate global PoD and add to dictionary
    global_pod_dict = get_global_pod(summary, seed)
    summary['global_pod_dict'] = global_pod_dict

    if no_compression:
        # Save as plain JSON without compression
        summary.to_json(path_to_summary, orient='index')
    else:
        # Save with zip compression
        summary.to_json(path_to_summary, orient='index', compression='zip')


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--fits-dir', type=str, help='list of probe .pkl files to process separated by spaces')
    parser.add_argument('--output', type=str, help='path to the output json')
    parser.add_argument('--seed', type=int, help='optional random seed for reproducibility')
    parser.add_argument('--no-compression', action='store_true', help='save output as plain JSON without compression')
    args = parser.parse_args()

    compress_output(args.fits_dir, args.output, args.seed, args.no_compression)


if __name__ == '__main__':
    main()
