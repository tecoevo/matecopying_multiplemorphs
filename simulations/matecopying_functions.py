import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import pandas as pd
from scipy.stats import binom
from scipy.stats import betabinom
from collections import Counter
from itertools import product
import matecopying_params as params

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

def generate_grid(n):
    grid = []
    for x in np.linspace(0, 1, n):
        for y in np.linspace(0, 1 - x, n - int(x * n)):
            z = 1 - x - y
            grid.append((x, y, z))
    return grid

def choose(m, array, probabilities):
    """ chooses an individual of one type out of 'm' types w.p. 'probabilities' and makes sure 
    that there is at least one individual of that type in 'array' """

    chosen = np.random.choice(range(m), p=probabilities)
    if array[chosen]<1:
        while array[chosen]<1:
            chosen = np.random.choice(range(m), p=probabilities)

    return chosen

def one_generation(m, b, c, y_t, q_values, n_population):
    """ simulating n_matings number of matings with a given copying probability"""

    q_proportion = q_values/sum(q_values)
    n_assort = np.round(q_proportion*n_population) # number of females that prefer each male morph (inherent prference)
    n_morphpop = np.round(y_t*n_population) # array of population of each morph
    
    for mating in range(int(n_matings)):
        # no copying
        if rnd.random()<(1-c):
            
            ### Pair choosing
            chosen_fem = choose(m,n_assort, q_proportion)
            chosen_male = choose(m,n_morphpop,y_t)
        
            if chosen_male!=chosen_fem:
                while(chosen_fem!=chosen_male):
                    chosen_fem = choose(m,n_assort, q_proportion)
                    chosen_male = choose(m,n_morphpop,y_t)

        
            # getting number of births
            q = q_values[chosen_male] # geting the male morph's quality
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births # adding births to the male morph population 
        
        else: # copying
            switching_probs = (y_t**b)/(y_t**b + (1-y_t)**b) # switching probabilities
            
            chosen_fem = choose(m,n_assort,q_proportion)
            # normalizing the probabilities
            switching_probs[chosen_fem] = 0
            prob_sum = sum(switching_probs)
            switching_probs[chosen_fem] = np.maximum(1-prob_sum,0)

            ### Pair choosing
            new_pref = np.random.choice(range(m), p=switching_probs)
            chosen_male = choose(m, n_morphpop, y_t)

            if chosen_male!=new_pref:
                while(new_pref!=chosen_male):
                    chosen_fem = choose(m,n_assort,q_proportion)
                    # normalizing the probabilities
                    switching_probs[chosen_fem] = 0
                    prob_sum = sum(switching_probs)
                    switching_probs[chosen_fem] = np.maximum(1-prob_sum,0)
            
                    new_pref = np.random.choice(range(m), p=switching_probs)
                    chosen_male = choose(m,n_morphpop,y_t)

            # getting number of births
            q = q_values[chosen_male] # geting the male morph's quality
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births # adding births to the male morph population 
        
        ### deaths
        
        n_deaths = n_births # keeping the population constant
        for death in range(n_deaths):
            dying_one = choose(m, n_morphpop, y_t)
            n_morphpop[dying_one] -= 1      
        

        n_population = sum(n_morphpop)
        y_t = n_morphpop/sum(n_morphpop)
        if(1 in y_t):
            break
    

    return y_t

def one_run(n_population, q_values, m, T, b, c, y_0):
    y_t = y_0
    y_hist = [y_0]
    for t in range(T):
        if(1 not in y_t):
            y_t = one_generation(m, b, c, y_t, q_values, n_population) # updating frequencies after one generation
            y_hist.append(y_t)
        else:
            y_hist.append(y_t)

    return y_hist

def simulation(n_population, n_runs, m, T, b, c, y_0, q_values, out_dir):
        
    all_runs = []
    for run in range(n_runs):
        y_hist = one_run(n_population, q_values, m, T, b, c, y_0)
        all_runs.append(y_hist)

    df = pd.DataFrame(all_runs)
    if m==3:
        df.to_pickle(f'{out_dir}/{float(y_0[0]):.2f}_{float(y_0[2]):.2f}_{float(c):.2f}.pkl')
    if m==2:
        df.to_pickle(f'{out_dir}/{float(y_0[1]):.2f}_{float(c):.2f}.pkl')
    
def run_simulations_2m(y,c, out_dir):
    simulation(n_population, n_runs, m, T, b, c, np.array([1-y,y]), q_values, out_dir)
    return f"Output written for y={y} and c={c}"

def run_simulations_3m(y,c, out_dir):
    simulation(n_population, n_runs, m, T, b, c, np.array(y), q_values, out_dir)
    return f"Output written for y={y} and c={c}"

def fix_fractions(x,y,c,input_dir):
    f = f'{float(x):.2f}_{float(y):.2f}_{float(c):.2f}.pkl'
    df = pd.read_pickle(f'{input_dir}/{f}')
    last_gen2 = np.sum([df[df.shape[1]-1][i][2] for i in range(df.shape[0])])/50
    last_gen0 = np.sum([df[df.shape[1]-1][i][0] for i in range(df.shape[0])])/50
    last_gen1 = np.sum([df[df.shape[1]-1][i][1] for i in range(df.shape[0])])/50
    return last_gen0, last_gen1, last_gen2