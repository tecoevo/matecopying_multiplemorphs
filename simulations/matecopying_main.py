import numpy as np
import pandas as pd
from itertools import product
from joblib import Parallel, delayed
import os
import itertools
import matecopying_params as params
from matecopying_functions import *

# Output directory
out_dir = "simulation_data/test/"
data_dir = f'{out_dir}/data'
out_file = f'{out_dir}/out.csv'
params_file = f'{out_dir}/params.csv'
comment_file = f'{out_dir}/comment.csv'

# Add any comment about this run of simulations, if needed
# Will be saved in a separate file in out_dir
comment = ""

### Import parameters and initial conditions
m = params.m
b = params.b
f = params.f 
copying_type = params.copying_type 
threshold_2m = params.threshold_2m 
threshold_3m = params.threshold_3m 
strong_anticonformity = params.strong_anticonformity 
n_population = params.n_population 
n_matings = params.n_matings 
n_runs = params.n_runs 
T = params.T 
y_range = params.y_range 
q_values = params.q_values 
c_range = params.c_range 


# Create directory if it doesn't exist
if not os.path.exists(data_dir):
    os.makedirs(out_dir)
    os.makedirs(data_dir)
else:
    print(f"Directory '{data_dir}' already exists.")


# Save params and comment to files in the output dir
params_dict = {
    "n_population": n_population,
    "n_matings": n_matings,
    "n_runs": n_runs,
    "m": m,
    "T": T,
    "b": b,
    "f": f,
    "threshold_2m": threshold_2m,
    "threshold_3m": threshold_3m,
    "strong_anticonformity": strong_anticonformity,
    "copying_type": copying_type,
    "y_range": y_range,
    "q_values": q_values,
    "c_range": c_range 
}

params_df = pd.DataFrame(list(params_dict.items()), columns=["Parameter", "Value"])
params_df.to_csv(params_file, index=False)
pd.DataFrame({"Comment": [comment]}).to_csv(comment_file, index=False)

### Run the simulations in parallel
if m==2:
    results = Parallel(n_jobs=5)(delayed(run_simulations_2m)(y,c, data_dir) for y,c in product(y_range, c_range))
    
    allfiles = os.listdir(data_dir)
    files = sorted(allfiles)
    params_lattice = list(itertools.product(y_range, c_range))
    data = []
    for f in files:
        coordinate = params_lattice[files.index(f)]
        df = pd.read_pickle(f'{data_dir}/{f}')

        # Average the outcome over all the runs
        x_coord = coordinate[1]
        y_coord = coordinate[0]
        z_value = np.mean([df[df.shape[1]-1][i][1] for i in range(df.shape[0])])
        data.append([x_coord, y_coord, z_value])

    df_save = pd.DataFrame(data, columns=["x", "y", "z"])
    df_save.to_csv(out_file, index=False)

if m==3:
    results = Parallel(n_jobs=5)(delayed(run_simulations_3m)(y,c, data_dir) for y,c in product(y_range, c_range))
    
    allfiles = os.listdir(data_dir)
    files = sorted(allfiles)
    x_coords = [point[0] for point in y_range]
    y_coords = [point[2] for point in y_range]

    for c in c_range:
        z_values =[fix_fractions(x,z,c,data_dir) for x,y,z in y_range]
        highest, middle, lowest = zip(*z_values)
        df_save = pd.DataFrame({'h':x_coords, 'l':y_coords, 'h_eq':highest, 'l_eq':lowest})
        df_save.to_csv(f'{out_dir}/out_c{c}.csv', index=False)