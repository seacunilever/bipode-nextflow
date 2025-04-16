#!/usr/bin/env python
import argparse
import pickle
import os
from pathlib import Path
from typing import Dict, Any, Union, List, Optional
import numpy as np
import pandas as pd


def process_data(input_file_path: Union[str, Path], path_to_output: Union[str, Path], testing_mode: bool = False) -> None:
    """
    Process raw count data into a stan-compatible format for each probe.

    Args:
        input_file_path: Path to pipeline input json file
        path_to_output: Path to the temporary directory where intermediate output will be stored
        testing_mode: If True, only process the first probe for testing purposes

    Returns:
        None
    """
    # Create directories for intermediate output within specified analysis directory
    output_path = Path(path_to_output)
    data_dir = output_path / "Data"
    fits_dir = output_path / "Fits"
    
    data_dir.mkdir(parents=True, exist_ok=True)
    fits_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    if not Path(input_file_path).exists():
        raise FileNotFoundError('json does not exist')
    
    df = pd.read_json(input_file_path, typ='series', orient='index')

    count_matrix = np.array(df['counts'], dtype='int')
    n_sample = count_matrix.shape[1]

    if 'total_count' not in df.index:
        total_count = np.sum(count_matrix, axis=0)
    else:
        total_count = np.array(df['total_count'], dtype='int')

    # Determine number of batches
    batch_index = np.array(df['batch_index'], dtype='int')
    n_batch = np.max(batch_index)
    n_treatment_batch = df['n_treatment_batch']

    # Calculate concentration index for each sample
    n = len(total_count)
    concentration = np.array(df['concentration'], dtype='float')
    unique_concentration = list(np.unique(concentration[concentration > 0]))
    n_conc = len(unique_concentration)
    concentration_index = np.zeros(n, dtype='int')
    for i, j in enumerate(concentration):
        if j in unique_concentration:
            concentration_index[i] = unique_concentration.index(j) + 1
    unique_concentration = np.log10(unique_concentration)

    for probe_index, (probe, probe_count) in enumerate(zip(df['probes'], count_matrix)):
        if testing_mode and probe_index > 0:
            break

        # Split counts by high and low
        low_count_index = np.where(probe_count <= 100)[0] + 1
        high_count_index = np.where(probe_count > 100)[0] + 1
        n_low_count = len(low_count_index)
        n_high_count = len(high_count_index)

        data: Dict[str, Any] = {
            'n_sample': n_sample,
            'n_treatment_batch': n_treatment_batch,
            'count': probe_count,
            'total_count': total_count,

            'n_batch': n_batch,
            'batch_index': batch_index,

            'n_conc': n_conc,
            'conc': unique_concentration,
            'conc_index': concentration_index,

            'n_low_count': n_low_count,
            'low_count_index': low_count_index,

            'n_high_count': n_high_count,
            'high_count_index': high_count_index,
        }

        with open(data_dir / f"{probe}.pkl", 'wb') as f:
            pickle.dump(data, f)


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, help='path to input data json')
    parser.add_argument('--analysis-dir', type=str, help='path to analysis directory')
    args = parser.parse_args()

    process_data(args.input_file, args.analysis_dir)


if __name__ == '__main__':
    main()