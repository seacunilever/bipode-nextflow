#!/usr/bin/env python
import argparse
import pickle
import os
import tarfile
import shutil
from pathlib import Path
from typing import Dict, Any, Union, List, Optional
import numpy as np
import pandas as pd
import glob


def create_manifest(batch_files: List[str], manifest_path: Path, tar_filename: str, batch_num: int) -> None:
    """Create a manifest file entry for a batch of files."""
    # If file doesn't exist, create it with header
    if not manifest_path.exists():
        with open(manifest_path, 'w') as f:
            f.write("batch\ttar_file\tprobes\n")

    # Get all probe names for this batch
    probe_names = [Path(file).stem for file in batch_files]
    # Join them with commas
    probes_str = ",".join(probe_names)

    with open(manifest_path, 'a') as f:
        f.write(f"{batch_num}\t{tar_filename}\t{probes_str}\n")


def create_tar_archive(files: List[str], output_path: Path) -> None:
    """Create a tar.gz archive from a list of files."""
    with tarfile.open(output_path, 'w:gz') as tar:
        for file in files:
            tar.add(file, arcname=Path(file).name)


def create_directory_archive(files: List[str], output_path: Path) -> None:
    """Create a directory containing the files."""
    output_path.mkdir(parents=True, exist_ok=True)
    for file in files:
        # Copy file to output directory
        shutil.copy2(file, output_path / Path(file).name)


def process_batches(data_dir: Path, prefix: str, batch_size: int, batch_mode: str, archive_mode: str = "tar") -> Path:
    """
    Process the pickle files into batches and create manifests and archives.

    Args:
        data_dir: Directory containing the pickle files
        prefix: Prefix for output files
        batch_size: Number of files per batch
        batch_mode: Either 'batch' or 'all' to control archiving behavior
        archive_mode: Either 'tar' or 'directory' to control how files are stored

    Returns:
        Path: Path to the created manifest file

    Raises:
        FileNotFoundError: If no pickle files are found in the data directory
        ValueError: If archive_mode is not 'tar' or 'directory'
    """
    if archive_mode not in ["tar", "directory"]:
        raise ValueError("archive_mode must be either 'tar' or 'directory'")

    # Get list of all pickle files
    files = sorted(glob.glob(str(data_dir / "*.pkl")))
    total_files = len(files)

    if not files:
        raise FileNotFoundError(f"No pickle files found in {data_dir}")

    # Create single manifest file
    manifest_path = Path(f"{prefix}.manifest.csv")
    # Clear the manifest file if it exists
    manifest_path.unlink(missing_ok=True)

    # Create single archive if batch_mode is 'all'
    if batch_mode == 'all':
        archive_name = f"{prefix}_batch0"
        if archive_mode == "tar":
            archive_path = Path(f"{archive_name}.tar.gz")
            create_tar_archive(files, archive_path)
        else:  # directory mode
            archive_path = Path(archive_name + "_dir")
            create_directory_archive(files, archive_path)

        # Process files in batches for manifest, even though we're using a single archive
        batch_num = 1
        for i in range(0, total_files, batch_size):
            batch_files = files[i:i + batch_size]
            create_manifest(batch_files, manifest_path, str(archive_path), batch_num)
            batch_num += 1
        return manifest_path

    # Process files in batches
    batch_num = 1
    for i in range(0, total_files, batch_size):
        batch_files = files[i:i + batch_size]
        batch_prefix = f"{prefix}_batch{batch_num}"

        if archive_mode == "tar":
            archive_path = Path(f"{batch_prefix}.tar.gz")
            create_tar_archive(batch_files, archive_path)
        else:  # directory mode
            archive_path = Path(batch_prefix)
            create_directory_archive(batch_files, archive_path)

        # Add entries to manifest
        create_manifest(batch_files, manifest_path, str(archive_path), batch_num)
        batch_num += 1

    return manifest_path


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
    parser.add_argument('--archive-mode', type=str, choices=['tar', 'directory'], default='tar',
                      help='archive mode: "tar" for tar.gz files, "directory" for directories')
    parser.add_argument('--prefix', type=str, help='prefix for output files')
    args = parser.parse_args()

    # Process the data into pickle files
    process_data(args.input_file, args.analysis_dir)

    # Process the pickle files into batches
    data_dir = Path(args.analysis_dir) / "Data"
    manifest_file = process_batches(data_dir, args.prefix, args.batch_size, args.batch_mode, args.archive_mode)


if __name__ == '__main__':
    main()
