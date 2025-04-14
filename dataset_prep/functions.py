import numpy as np
import pandas as pd
import yaml


def load_yaml_file(file_path) -> dict:
    """
    Opens file and returns contents as string.

    Accepts:
        file_path (str) - string obtained from upload button

    Returns:
        yaml_dict (dict) - string to be displayed in mapper editor box
    """

    with open(file_path) as stream:
        try:
            yaml_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return yaml_dict


def convert_meta_data(meta: pd.DataFrame, meta_mapper_dict: dict) -> pd.DataFrame:
    """
    Builds a pandas DataFrame containing meta data for internal use
    in the app.

    Accepts:
        meta (pd.DataFrame) - pandas DataFrame of meta data
        meta_data_dict (dict) - dictionary containing possible meta data key names for conversion

    Returns:
        df  (pd.DataFrame) - pandas DataFrame of meta data
    """

    # Populate columns from meta data and mapper dict
    df = pd.DataFrame()
    for key in meta_mapper_dict:
        for i in meta_mapper_dict[key]:

            if isinstance(i, str):
                if i in meta.columns:
                    df[key] = meta[i]
                    break

            elif isinstance(i, dict):
                subkey = list(i.keys())[0]
                if subkey in meta.columns:
                    # Check if MetaDataConversionUtilities class has a method corresponding
                    # to the required function named defined in the meta data mapper yaml
                    if i[subkey] in dir(MetaDataConversionUtilities):
                        func = getattr(MetaDataConversionUtilities, i[subkey])
                        df[key] = func(meta[subkey], meta)
                        break
                    else:
                        raise AttributeError(f'No function defined for {i[subkey]}')

            else:
                raise TypeError(f'No logic defined to handle key of type {type(i)}')

    # Remove any rows with nans in the Cell type column
    if 'Cell type' in df.columns:
        df = df[~df['Cell type'].isna()]

    return df


def filter_percent_mapped_reads(df: pd.DataFrame, minimum_percent_mapped_reads: (int, float)) -> pd.DataFrame:
    """
    Filters out samples below the specified minimum percentage of mapped reads.

    Accepts:
        df (pd.DataFrame) - pandas DataFrame of sample meta data
        minimum_percent_mapped_reads (int, float) - threshold to filter percent mapped reads

    Returns:
        df (pd.DataFrame) - filtered meta data
    """

    df = df[df['Percent mapped reads'] >= minimum_percent_mapped_reads]

    return df


def filter_total_mapped_reads(df: pd.DataFrame, minimum_total_mapped_reads: (int, float)) -> pd.DataFrame:
    """
    Filters out samples below the specified minimum total mapped reads.

    Accepts:
        df (pd.DataFrame) - pandas DataFrame of sample meta data
        minimum_total_mapped_reads (int, float) - threshold to filter total mapped reads

    Returns:
        df (pd.DataFrame) - filtered meta data
    """

    df = df[df['Num. mapped reads'] >= minimum_total_mapped_reads]

    return df


def write_bifrost_input(meta, filter_dict, counts_table, config_dict, output_directory):
    """
    Applies filters to DataFrame and returns BIFROST HTTr pipeline input.

    Accepts:
        meta (pd.DataFrame) - DataFrame of meta data
        filters (dict) - dictionary of filters to apply
        counts_table (pd.DataFrame) - DataFrame of meta data
        config_dict (dict) - dictionary of configuration settings
        file_path (str) - path to location of saved input file

    Returns:
        None
    """

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
    """
    Generates BIFROST inputs from the provided meta DataFrame, counts DataFrame
    and config dict.

    Accepts:
        meta (pd.DataFrame) - DataFrame of meta data
        counts_table (pd.DataFrame) - DataFrame of meta data
        config_dict (dict) - dictionary of configuration settings
        output_directory (str) - directory into which files will be written

    Returns:
        None
    """

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