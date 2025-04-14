import os
import shutil
import pickle

from bin.split_data import process_data
from bin.run_conc_response_analysis import run_concentration_response_analysis, fit_model, gen_plotting_data
from bin.compress_output import compress_output

input_file = f'dataset_prep/bifrost_inputs/BIFROST_input_Nitrofurantoin_HepG2.json'
path_to_stan_code = f'model/BIFROST_HTTr_beta_logistic_batch.stan'
temp_dir = f'temp'


if __name__ == '__main__':

    # Remove temp directory, and reinstall
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
    os.mkdir(f'{temp_dir}/Data')
    os.mkdir(f'{temp_dir}/Fits')

    # Attempt to process input file
    process_data(input_file, temp_dir)

    # # Load model and first probe-specific input
    probes = os.listdir(f'{temp_dir}/Data')
    data_paths = [f'{temp_dir}/Data/{i}' for i in probes]
    run_concentration_response_analysis(data_paths, path_to_stan_code, 4, fit_dir=temp_dir)

    # Attempt to compress output
    compress_output(f'{temp_dir}/Fits', 'temp/summary.json.zip')
