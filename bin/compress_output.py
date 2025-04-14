#!/usr/bin/env python

import argparse
import os
import numpy as np
import pandas as pd
import pickle


def get_confidence_threshold_probability_density(x: np.ndarray):
    """
    Evaluates the probability density for the function for
    defined to describe uncertainty in CDS threshold.

    Accepts:
        x (np.ndarray): array of values at which to calculate density

    Returns:
        dq (np.ndarray) - corresponding probability density
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


def get_minimum_pod_means(pod_means: np.ndarray, cds: np.ndarray, cds_thresholds: np.ndarray, max_conc: float):
    """
    Computes minimum PoD means for each value in the specified array

    Accepts:
        stats (pd.Series) - dictionary-like structure with PoD means and CDS
        cds_thresholds (np.ndarray) - thresholds on which to filter probes

    Returns:
        min_means (np.ndarray) - array of minimum mean values
        min_probes (np.ndarray) - array of corresponding probe IDs
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


def get_global_pod(df: pd.Series):
    """
    Calculates and returns global PoD.

    Accepts:
        df (pd.Series): BIFROST summary

    Returns:
        results (dict): dictionary of global PoD-related stats
    """

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


def compress_output(analysis_dir, path_to_summary):
    """
    Compresses intermediate output into a single pandas Dataframe.

    Accepts:
        analysis_dir - name of parent directory for intermediate output
        path_to_summary - path to summary file

    Returns:
        None
    """

    # Determine probe IDs
    data_files = os.listdir(f'{analysis_dir}')
    probes = np.array([os.path.splitext(file)[0] for file in data_files])

    # Create empty pandas series
    summary = pd.Series(dtype='object')

    # Extract details inputs universal to all chemicals/probes
    data = pickle.load(open(f'{analysis_dir}/{probes[0]}.pkl', 'rb'))
    for key in ['n_samp',
                'n_sample', 'n_treatment_batch', 'total_count', 'n_batch', 'batch_index',
                'n_conc', 'conc', 'conc_index', 'max_conc',
                ]:
        summary[key] = data[key]

    summary['probes'] = probes

    # Extract probe-specific information
    for probe in probes:
        data = pickle.load(open(f'{analysis_dir}/{probe}.pkl', 'rb'))

        summary[probe] = pd.Series(dtype='object')
        summary[probe]['diagnostics'] = data['diagnostics']
        summary[probe]['parameters'] = data['parameters']

        for par in data['fit'].index:
            summary[probe][par] = data['fit'][par]
        summary[probe]['count'] = data['count']

    # Calculate global PoD and add to dictionary
    global_pod_dict = get_global_pod(summary)
    summary['global_pod_dict'] = global_pod_dict

    summary.to_json(path_to_summary, orient='index', compression='zip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fits-dir', type=str, help='list of probe .pkl files to process separated by spaces')
    parser.add_argument('--output', type=str, help='path to the output json')
    args = parser.parse_args()

    compress_output(args.fits_dir, args.output)