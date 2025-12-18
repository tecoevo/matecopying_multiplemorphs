import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import pandas as pd
from scipy.stats import binom
from scipy.stats import betabinom
from collections import Counter
from itertools import product
from joblib import Parallel, delayed
import os
import itertools

#script_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(script_dir)

import matecopying_params as params
from matecopying_functions import *

### import parameters
n_population = params.n_population # male population, total population is twice this number
n_matings = params.n_matings # number of matings in one generation
n_runs = params.n_runs # number of independent runs
m = params.m # number of male morphs
T = params.T # number of generations
b = params.b # extent of conformism (beta)
y_range = params.y_range # range of initial frequencies
q_values = params.q_values # male morph qualities
c_range = params.c_range # range of copying probabilities
copying_type = params.copying_type
factor = params.factor
threshold_2m = params.threshold_2m # threshold for copying type 3, m=2
threshold_3m = params.threshold_3m # threshold for copying type 3, m>2
BR_type = params.BR_type # type of anticonformity function

# directory to save simulation data
out_dir = "simulation_data/2m_real_f2/"
data_dir = f'{out_dir}/data'
out_file = f'{out_dir}/out.csv'
params_file = f'{out_dir}/params.csv'
comment_file = f'{out_dir}/comment.csv'
comment = ""
#, thresholds = (0.7, 0.5)

# create directory if it doesn't exist
if not os.path.exists(data_dir):
    os.makedirs(out_dir)
    os.makedirs(data_dir)
else:
    print(f"Directory '{data_dir}' already exists.")


# save params and comment to files in the output dir
params_dict = {
    "n_population": n_population,
    "n_matings": n_matings,
    "n_runs": n_runs,
    "m": m,
    "T": T,
    "b": b,
    "y_range": y_range,
    "q_values": q_values,
    "c_range": c_range,
    "copying_type": copying_type,
    "factor": factor,
    "threshold_2m": threshold_2m,
    "threshold_3m": threshold_3m,
    "BR_type": BR_type
}

params_df = pd.DataFrame(list(params_dict.items()), columns=["Parameter", "Value"])
params_df.to_csv(params_file, index=False)

pd.DataFrame({"Comment": [comment]}).to_csv(comment_file, index=False)

### run the simulations in parallel
if m==2:
    results = Parallel(n_jobs=5)(delayed(run_simulations_2m)(y,c, data_dir) for y,c in product(y_range, c_range))
    #print(results)
    allfiles = os.listdir(data_dir)
    files = sorted(allfiles)
    params_lattice = list(itertools.product(y_range, c_range))
    data = []
    for f in files:
        coordinate = params_lattice[files.index(f)]
        df = pd.read_pickle(f'{data_dir}/{f}')

        # average the outcome over all the runs
        x_coord = coordinate[1]
        y_coord = coordinate[0]
        z_value = np.mean([df[df.shape[1]-1][i][1] for i in range(df.shape[0])])
        data.append([x_coord, y_coord, z_value])

    df_save = pd.DataFrame(data, columns=["x", "y", "z"])
    df_save.to_csv(out_file, index=False)

if m==3:
    results = Parallel(n_jobs=5)(delayed(run_simulations_3m)(y,c, data_dir) for y,c in product(y_range, c_range))
    #print(results)
    allfiles = os.listdir(data_dir)
    files = sorted(allfiles)
    x_coords = [point[0] for point in y_range]
    y_coords = [point[2] for point in y_range]

    for c in c_range:
        z_values =[fix_fractions(x,z,c,data_dir) for x,y,z in y_range]
        highest, middle, lowest = zip(*z_values)
        df_save = pd.DataFrame({'h':x_coords, 'l':y_coords, 'h_eq':highest, 'l_eq':lowest})
        df_save.to_csv(f'{out_dir}/out_c{c}.csv', index=False)