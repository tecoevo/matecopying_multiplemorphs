import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import pandas as pd
from scipy.stats import binom
from scipy.stats import betabinom
import matplotlib
from collections import Counter
import sys
import argparse

# read initial values from input file (first argument)
parser = argparse.ArgumentParser()
parser.add_argument('--param', type=float, default=0.0, nargs="?")
args = parser.parse_args()

print("START", file = sys.stderr) 
print(sys.argv, file= sys.stderr)

# parameter values
n_population = 100 # male population, total population is twice this number
n_matings = n_population # number of matings in one generation
n_runs = 100 # number of runs
n_obs = 10 # number of observations
m = 2 # number of male morphs
T = 100 # number of generations
b = 20 # extent of conformism
u = n_obs/n_population # training period
search = 2 # search parameter

# range of gamma values
c_range = np.around(np.arange(0,1.001,0.01), 2)

def choose(m, array, probabilities):
    """ chooses an individual of one type out of 'm' types w.p. 'probabilities' and makes sure that there is at least one individual of that type in 'array' """

    chosen = np.random.choice(range(m), p=probabilities)
    if array[chosen]<1:
        while array[chosen]<1:
            chosen = np.random.choice(range(m), p=probabilities)

    return chosen

def one_generation(n_obs, m, b, c, y_t, q_values, n_population):
    """ one generation has a training period with n_obs number of matings, in which there is no copying behaviour. The rest of the matings have copying behaviour"""

    q_proportion = q_values/sum(q_values)
    n_assort = np.round_(q_proportion*n_population) # number of females that prefer each male morph (inherent prference)
    n_morphpop = np.round_(y_t*n_population) # array of population of each morph
    
    
    ###########################################################################
    ### training period
    ###########################################################################
    
    training = [] # array of male morphs chosen for each mating
    
    for obs in range(n_obs):
                
        ### births
        # the female searches for her preference 'search' number of times
        
        female_pref = choose(m, n_assort, q_proportion) # choose female type to mate

        chosen_one = choose(m, n_morphpop, y_t) # choose a male morph
        for s in range(search): # search for the right male morph
            
            if chosen_one == female_pref: 
                break
            else:
                chosen_one = choose(m, n_morphpop, y_t) # choose another male morph

        # getting number of births
        q = q_values[chosen_one] # geting the male morph's quality
        if q==1:
            n_births = 1
        else:
            alpha = (q-1)/(6-q) # parameter for choosing number of offspring
            n_births = betabinom.rvs(5, alpha, 1) # sampling from a distribution to get number of children
            
        n_morphpop[chosen_one] += n_births # adding births to the male morph population 
        
        ### deaths
        
        n_deaths = n_births # keeping the population constant
        for death in range(n_deaths):
            dying_one = choose(m, n_morphpop, y_t)
            n_morphpop[dying_one] -= 1
        
        training.append(chosen_one) # add choice to the training dataset

    ###########################################################################
    ### copying period
    ###########################################################################
        
    for mating in range(int(n_population - n_obs)):
        
        # fraction of male morph appearance in training set
        count = Counter(training)
        p = np.zeros(m)
        for i in range(m):
            p[i] = count[i]/n_obs

        ### births
        # initial female morph distribution: number of females is equal to the number of males. The frequency distribution is according to quality values, i.e. more females prefer the higher quality male
        
        # no copying
        if rnd.random()<(1-c):
            female_pref = choose(m, n_assort, q_proportion) # choose female type to mate

            chosen_one = choose(m, n_morphpop, y_t) # choose a male morph
            for s in range(search): # search for the right male morph
                
                if chosen_one == female_pref: 
                    break
                else:
                    chosen_one = choose(m, n_morphpop, y_t) # choose another male morph

            # getting number of births
            q = q_values[chosen_one] # geting the male morph's quality
            if q==1:
                n_births = 1
            else:
                alpha = (q-1)/(6-q) # parameter for choosing number of offspring
                n_births = betabinom.rvs(5, alpha, 1) # sampling from a distribution to get number of children
            
            n_morphpop[chosen_one] += n_births # adding births to the male morph population 
        
        else: # copying
            switching_probs = 1/(1+np.exp(-b*(p-0.5))) # switching probabilities

            # normalizing the probabilities
            if sum(switching_probs)<1:
                switching_probs[female_pref] = 0
                prob_sum = sum(switching_probs)
                switching_probs[female_pref] = 1-sum(switching_probs)
            else:
                switching_probs = switching_probs/sum(switching_probs)

            female_pref = choose(m, n_assort, switching_probs) # choose female type to mate
            chosen_one = choose(m, n_morphpop, y_t) # choose a male morph
            for s in range(search):
                if chosen_one==female_pref:
                    break
                else:
                    chosen_one = choose(m, n_morphpop, y_t)
              

            # getting number of births
            q = q_values[chosen_one] # geting the male morph's quality
            if q==1:
                n_births = 1
            else:
                alpha = (q-1)/(6-q) # parameter for choosing number of offspring
                n_births = betabinom.rvs(5, alpha, 1) # sampling from a distribution to get number of children
            
            n_morphpop[chosen_one] += n_births # adding births to the male morph population 
        
            # updating female population preferences           
            if chosen_one != female_pref: # if switching happens
                n_assort[chosen_one] += 1
                n_assort[female_pref] -= 1
                    
        ### deaths
        
        n_deaths = n_births # keeping the population constant
        for death in range(n_deaths):
            dying_one = choose(m, n_morphpop, y_t)
            n_morphpop[dying_one] -= 1      
        

        n_population = sum(n_morphpop)
        
    y_t = n_morphpop/sum(n_morphpop)

    return y_t


def one_run(n_population, n_obs, q_values, m, T, b, c, y_0):
    y_t = y_0
    y_hist = [y_0]
    for t in range(T):
        y_t = one_generation(n_obs, m, b, c, y_t, q_values, n_population) # updating frequencies after one generation
        y_hist.append(y_t)

    return y_hist


def simulation(n_population, n_obs, n_runs, m, T, b, c, y_0, q_values):
        
    all_runs = []
    for run in range(n_runs):
        y_hist = one_run(n_population, n_obs, q_values, m, T, b, c, y_0)
        all_runs.append(y_hist)

    df = pd.DataFrame(all_runs)
    df.to_pickle('simulation_data_{y}_{c}.pkl'.format(y="{:.2f}".format(y_0[1]), c="{:.2f}".format(c)))
      
y = float(args.param)
y_0 = np.array([1-y,y ])    
q_values = np.array([2.3,2])

for c in c_range:
    print('gamma = ', c)
    print('y_0 for lower quality male = ',y)
    simulation(n_population, n_obs, n_runs, m, T, b, c, y_0, q_values)

# saving params file
params = np.array(['population', 'observations', 'runs', 'morphs', 'generations', 'quality values'])
param_vals = np.array([n_population, n_obs, n_runs, m, T, q_values])
params_df = pd.DataFrame([params, param_vals])
params_df.to_csv('parameters.csv')

