import os
import numpy as np
import pandas as pd

from functions import load_yaml_file, convert_meta_data, generate_bifrost_inputs

# Load meta data file and convert in to common format
meta_raw = pd.read_csv('data/Example_Meta_Data.csv')
meta_data_mapper = load_yaml_file('data/sers_meta_data_mapper.yml')
meta = convert_meta_data(meta_raw, meta_data_mapper)

# Load counts table
counts = pd.read_csv('data/Example_Counts.csv')

# Load configuration file specifying processing choices
config_dict = load_yaml_file('data/example_config_file.yaml')

# Make directory for storing outputs
output_dir = 'bifrost_inputs'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
generate_bifrost_inputs(meta, counts, config_dict, output_dir)

# Reduce HepG2 datasets to 5 probes only (for testing pipeline on small datasets)
for i in os.listdir(output_dir):
    if 'HepG2' in i:
        df = pd.read_json(f'{output_dir}/{i}', orient='index', typ='series')
        print(df['probes'][:5])
        index = np.random.choice(len(df['probes']), size=5, replace=False)
        df['probes'] = np.array(df['probes'])[index]
        df['counts'] = np.array(df['counts'], dtype='int')[index]
        df.to_json(f'{output_dir}/{i}', orient='index')
