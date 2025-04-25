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


# directory to save simulation data
out_dir = "simulation_data/final_test_2m/"
out_file = "simulation_data/final_test_2m/out.csv"

# create directory if it doesn't exist
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
else:
    print(f"Directory '{out_dir}' already exists.")

### run the simulations in parallel
if m==2:
    results = Parallel(n_jobs=5)(delayed(run_simulations_2m)(y,c) for y,c in product(y_range, c_range))
    #print(results)
    allfiles = os.listdir(out_dir)
    files = sorted(allfiles)
    params_lattice = list(itertools.product(y_range, c_range))
    data = []
    for f in files:
        coordinate = params_lattice[files.index(f)]
        df = pd.read_pickle(f'{out_dir}/{f}')

        # average the outcome over all the runs
        x_coord = coordinate[1]
        y_coord = coordinate[0]
        z_value = np.mean([df[df.shape[1]-1][i][1] for i in range(df.shape[0])])
        data.append([x_coord, y_coord, z_value])

    df_save = pd.DataFrame(data, columns=["x", "y", "z"])
    df_save.to_csv(out_file, index=False)

if m==3:
    results = Parallel(n_jobs=5)(delayed(run_simulations_3m)(y,c) for y,c in product(y_range, c_range))
    #print(results)
    allfiles = os.listdir(out_dir)
    files = sorted(allfiles)
    x_coords = [point[0] for point in y_range]
    y_coords = [point[2] for point in y_range]

    for c in c_range:
        z_values =[fix_fractions(x,z,c,out_dir) for x,y,z in y_range]
        highest, middle, lowest = zip(*z_values)
        df_save = pd.DataFrame({'h':x_coords, 'l':y_coords, 'h_eq':highest, 'l_eq':lowest})
        df_save.to_csv(f'{out_dir}/out_c{c}.csv', index=False)