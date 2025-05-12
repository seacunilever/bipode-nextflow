#!/usr/bin/env python
import argparse
import pickle
import os
import tarfile
from pathlib import Path
from typing import Dict, Any, Union, List, Optional
import numpy as np
import pandas as pd
import glob


def create_manifest(batch_files: List[str], manifest_path: Path) -> None:
    """Create a manifest file for a batch of files."""
    with open(manifest_path, 'w') as f:
        for file in batch_files:
            # Remove .pkl extension as in original script
            f.write(f"{Path(file).stem}\n")


def create_tar_archive(files: List[str], output_path: Path) -> None:
    """Create a tar.gz archive from a list of files."""
    with tarfile.open(output_path, 'w:gz') as tar:
        for file in files:
            tar.add(file, arcname=Path(file).name)


def process_batches(data_dir: Path, prefix: str, batch_size: int, batch_mode: str) -> None:
    """
    Process the pickle files into batches and create manifests and archives.

    Args:
        data_dir: Directory containing the pickle files
        prefix: Prefix for output files
        batch_size: Number of files per batch
        batch_mode: Either 'batch' or 'all' to control archiving behavior
    """
    # Get list of all pickle files
    files = sorted(glob.glob(str(data_dir / "*.pkl")))
    total_files = len(files)

    if not files:
        return

    # Create single tar.gz if batch_mode is 'all'
    if batch_mode == 'all':
        create_tar_archive(files, Path(f"{prefix}_batch0.tar.gz"))

    # Process files in batches
    batch_num = 1
    for i in range(0, total_files, batch_size):
        batch_files = files[i:i + batch_size]
        batch_prefix = f"{prefix}_batch{batch_num}"

        # Create manifest file for this batch
        manifest_path = Path(f"{batch_prefix}.manifest.csv")
        create_manifest(batch_files, manifest_path)

        # Create individual tar file if in batch mode
        if batch_mode == 'batch':
            create_tar_archive(batch_files, Path(f"{batch_prefix}.tar.gz"))

        batch_num += 1


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
    parser.add_argument('--batch-size', type=int, default=0, help='number of files per batch')
    parser.add_argument('--batch-mode', type=str, choices=['batch', 'all'], default='all',
                      help='batch mode: "batch" for individual archives, "all" for single archive')
    parser.add_argument('--prefix', type=str, help='prefix for output files')
    args = parser.parse_args()

    # Process the data into pickle files
    process_data(args.input_file, args.analysis_dir)

    # Process the pickle files into batches
    data_dir = Path(args.analysis_dir) / "Data"
    process_batches(data_dir, args.prefix, args.batch_size, args.batch_mode)


if __name__ == '__main__':
    main()
