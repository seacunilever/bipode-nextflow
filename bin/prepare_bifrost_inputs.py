#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import yaml
import argparse
from pathlib import Path
import itertools

def load_yaml_file(file_path) -> dict:
    """Opens file and returns contents as string."""
    with open(file_path) as stream:
        try:
            yaml_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return yaml_dict

def convert_meta_data(meta: pd.DataFrame, meta_mapper_dict: dict) -> pd.DataFrame:
    """Builds a pandas DataFrame containing meta data for internal use."""
    # Populate columns from meta data and mapper dict
    df = pd.DataFrame()
    for key in meta_mapper_dict:
        for i in meta_mapper_dict[key]:
            if isinstance(i, str):
                if i in meta.columns:
                    df[key] = meta[i]
                    break
            else:
                raise TypeError(f'No logic defined to handle key of type {type(i)}')

    # Remove any rows with nans in the Cell type column
    if 'Cell type' in df.columns:
        df = df[~df['Cell type'].isna()]

    return df

def filter_percent_mapped_reads(df: pd.DataFrame, minimum_percent_mapped_reads: (int, float)) -> pd.DataFrame:
    """Filters out samples below the specified minimum percentage of mapped reads."""
    df = df[df['Percent mapped reads'] >= minimum_percent_mapped_reads]
    return df

def filter_total_mapped_reads(df: pd.DataFrame, minimum_total_mapped_reads: (int, float)) -> pd.DataFrame:
    """Filters out samples below the specified minimum total mapped reads."""
    df = df[df['Num. mapped reads'] >= minimum_total_mapped_reads]
    return df

def write_bifrost_input(meta, filter_dict, counts_table, config_dict, output_directory):
    """Applies filters to DataFrame and returns BIFROST HTTr pipeline input."""
    # Filter meta data
    test_substance_mask = meta['Test substance'] == filter_dict['Test substance']
    control_mask = meta['Concentration'] == 0
    mask = (test_substance_mask ^ control_mask)
    for key in filter_dict:
        if key not in ['Test substance', 'N/A']:
            additional_mask = meta[key] == filter_dict[key]
            mask = mask & additional_mask
    df = meta[mask]

    # Apply global filters
    df = filter_percent_mapped_reads(df, config_dict['Minimum percent mapped reads'])
    df = filter_total_mapped_reads(df, config_dict['Minimum number mapped reads'])

    # Apply specific filters
    for key in config_dict['Specific filters']:
        for value in config_dict['Specific filters'][key]:
            df = df[df[key] != value]

    # Define Stan variables
    concentration = df['Concentration'].values
    treatment_mask = concentration > 0
    unique_batches = list(df[config_dict['Batch key']].unique())
    batch_index = np.array([unique_batches.index(k) + 1 for k in df[config_dict['Batch key']]], dtype='int')
    n_treatment_batch = len(np.unique(batch_index[treatment_mask]))

    # Extract counts matrix
    probes = counts_table[counts_table.columns[0]].values
    counts = counts_table[df['Sample ID']].values.astype('int')

    # Filter probes if median or mean raw count is below
    treatment_mask = concentration > 0
    probe_to_retain_mask = np.array([
        True if (np.mean(k[treatment_mask]) > config_dict['Minimum average treatment count'] and
                 np.median(k[treatment_mask]) > config_dict['Minimum average treatment count'])
        else False for k in counts])
    probes = probes[probe_to_retain_mask]
    counts = counts[probe_to_retain_mask]

    # Write BIFROST input as dictionary and write to file
    bifrost_input = pd.Series({
        'test_substance': filter_dict['Test substance'],
        'cell_type': filter_dict['Cell type'],
        'probes': probes,
        'counts': counts,
        'batch_index': batch_index,
        'concentration': concentration,
        'n_treatment_batch': n_treatment_batch,
    })

    s = "".join(ch for ch in filter_dict['Test substance'] if ch.isalnum())
    for key in filter_dict:
        if key not in ['Test substance', 'N/A']:
            s += f'_{"".join(ch for ch in filter_dict[key] if ch.isalnum())}'
    file_path = f'{output_directory}/BIFROST_input_{s}.json'

    bifrost_input.to_json(file_path, orient='index')

def generate_bifrost_inputs(meta: pd.DataFrame, counts_table: pd.DataFrame, config_dict: dict, output_directory):
    """Generates BIFROST inputs from the provided meta DataFrame, counts DataFrame and config dict."""
    test_substances = config_dict['Test substances']
    cell_types = config_dict['Cell types']

    for i, test_substance in enumerate(test_substances):
        if config_dict['Additional divider'] == 'N/A':
            for cell_type in cell_types:
                filter_dict = {'Test substance': test_substance, 'Cell type': cell_type}
                write_bifrost_input(meta, filter_dict, counts_table, config_dict, output_directory)
        else:
            pairs = itertools.product(cell_types, meta[config_dict['Additional divider']].unique())
            for pair in pairs:
                filter_dict = {'Test substance': test_substance,
                               'Cell type': pair[0],
                               config_dict['Additional divider']: pair[1]}
                write_bifrost_input(meta, filter_dict, counts_table, config_dict, output_directory)

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare Bifrost inputs from meta data and counts')
    parser.add_argument('--meta-data', required=True, help='Path to meta data CSV file')
    parser.add_argument('--meta-mapper', required=True, help='Path to meta data mapper YAML file')
    parser.add_argument('--counts', required=True, help='Path to counts CSV file')
    parser.add_argument('--substances-cell-types', required=True, help='Path to substances and cell types YAML file')
    parser.add_argument('--additional-divider', default='N/A', help='Additional field to use for dividing data')
    parser.add_argument('--batch-key', default='Exposure plate ID', help='Field to use as batch key in the BIFROST model')
    parser.add_argument('--min-percent-mapped-reads', type=float, default=50.0, help='Minimum percentage of mapped reads required')
    parser.add_argument('--min-num-mapped-reads', type=int, default=100000, help='Minimum number of mapped reads required')
    parser.add_argument('--min-avg-treatment-count', type=float, default=5.0, help='Minimum average treatment count required')
    parser.add_argument('--specific-filters', default='{}', help='Additional specific filters to apply (JSON string)')
    parser.add_argument('--output-dir', default='bifrost_inputs', help='Directory to store outputs')
    parser.add_argument('--test-probes', type=int, default=5, help='Number of probes to sample for testing')
    parser.add_argument('--test-pattern', default='HepG2', help='Pattern to match files for testing')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load meta data file and convert in to common format
    meta_raw = pd.read_csv(args.meta_data)
    meta_data_mapper = load_yaml_file(args.meta_mapper)
    meta = convert_meta_data(meta_raw, meta_data_mapper)

    # Load counts table
    counts = pd.read_csv(args.counts)

    # Load substances and cell types
    substances_cell_types = load_yaml_file(args.substances_cell_types)

    # Create config dictionary from arguments
    config_dict = {
        'Test substances': substances_cell_types['Test substances'],
        'Cell types': substances_cell_types['Cell types'],
        'Additional divider': args.additional_divider,
        'Batch key': args.batch_key,
        'Minimum percent mapped reads': args.min_percent_mapped_reads,
        'Minimum number mapped reads': args.min_num_mapped_reads,
        'Minimum average treatment count': args.min_avg_treatment_count,
        'Specific filters': eval(args.specific_filters)  # Convert JSON string to dict
    }

    # Make directory for storing outputs
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    generate_bifrost_inputs(meta, counts, config_dict, output_dir)

    # Reduce datasets to specified number of probes (for testing pipeline on small datasets)
    for i in os.listdir(output_dir):
        if args.test_pattern in i:
            df = pd.read_json(f'{output_dir}/{i}', orient='index', typ='series')
            print(df['probes'][:args.test_probes])
            index = np.random.choice(len(df['probes']), size=args.test_probes, replace=False)
            df['probes'] = np.array(df['probes'])[index]
            df['counts'] = np.array(df['counts'], dtype='int')[index]
            df.to_json(f'{output_dir}/{i}', orient='index')

if __name__ == '__main__':
    main()
