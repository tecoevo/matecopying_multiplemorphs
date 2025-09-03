import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import pandas as pd
from scipy.stats import binom
from scipy.stats import betabinom
from collections import Counter
from itertools import product
import params as params

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
copying_type = params.copying_type # type of mate copying to determine copying probabilities
factor = params.factor # modulates the extent of conformity/anticonfrmity, higher means less conformity
threshold_2m = params.threshold_2m # threshold for copying type 3, m=2
threshold_3m = params.threshold_3m # threshold for copying type 3, m>2

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

def compute_D(greater_vals, b, factor):
    #if b<0: 
    #    print("Invalid value of b, assuming b=1")
    #    b=1
    if b>=1:
        y_max = np.max(greater_vals)
        D_max = (1/y_max - 1)*sum(greater_vals)
        return D_max/factor
    elif b<1:
        D_min = -sum(greater_vals)
        return D_min/factor
    
def copying_probabilities(y_t, b, factor, copying_type=1):
    """ computes switching probabilities based on the type of mate copying: 
    1: continuous copying function
    2: discrete copying function 
    3: custom copying function (more realistic)
    Output: array probabilities of switching to a given morph (may or may not
    depend on the chosen female)"""

    if copying_type==1:
        if b>=1: # not normalised, but will be normalised later for this case
            probs = (y_t**b)/(y_t**b + (1-y_t)**b)
        else: # normalised
            if b < 0:
                probs = np.where((y_t == 0) | (y_t == 1), y_t, (y_t**b) / (y_t**b + (1-y_t)**b)) 
            else:
                probs = (y_t**b)/(y_t**b + (1-y_t)**b)
            
            # normalizing the probabilities
            prob_sum = sum(probs)
            probs = probs/prob_sum
            
    elif copying_type==2:
        if b==1: probs = y_t # neutral copying
        else: 
            majority = 1/len(y_t)
            vec = np.array(y_t)
            
            # partition the frequencies
            greater_idx = np.where(vec > majority)[0]
            equal_idx = np.where(vec == majority)[0]
            less_idx = np.where(vec < majority)[0]

            greater_vals = vec[greater_idx]
            equal_vals = vec[equal_idx]
            less_vals = vec[less_idx]

            # compute value of D
            if len(greater_vals)==0: D = 0
            else: D = compute_D(greater_vals, b, factor)

            # convert to copying probabilities
            # greater than majority
            greater_vals = greater_vals + (greater_vals/sum(greater_vals))*D
            # less than majority
            inv_y = np.zeros_like(less_vals, dtype=float)
            mask = np.abs(less_vals) > 1e-12
            inv_y[mask] = 1 / less_vals[mask]

            if np.sum(inv_y) == 0:
                factors = np.zeros_like(inv_y)
            else:
                factors = inv_y / np.sum(inv_y)
                
            less_vals = less_vals - factors*D
            #inv_y = 1/less_vals
            #factors = inv_y / np.sum(inv_y)
            #less_vals = less_vals - factors*D
            
            # construct vector of copying probabilities
            probs = np.empty_like(vec)
            probs[greater_idx] = greater_vals
            probs[equal_idx] = equal_vals
            probs[less_idx] = less_vals
            probs = np.where(probs < 1e-12, 0, probs)

            #normalize
            prob_sum = sum(probs)
            probs = probs/prob_sum

    elif copying_type==3:
        majority = 1/m
        vec = np.array(y_t)
    
        # partition the frequencies 
        greater_idx = np.where(vec > majority)[0]
        equal_idx = np.where(vec == majority)[0]
        less_idx = np.where(vec < majority)[0]

        greater_vals = vec[greater_idx]
        equal_vals = vec[equal_idx]
        less_vals = vec[less_idx]

        # compute value of D
        if len(greater_vals)==0: 
            D_conf = 0
            D_anti = 0
        else: 
            D_conf = compute_D(greater_vals, b=2, factor = factor)
            D_anti = compute_D(greater_vals, b=-2, factor = factor)

        ##### testing
        #print([D_conf, D_anti, y_t])

        if m==2:
            threshold = threshold_2m
        else: # for more than 2 morphs
            threshold = threshold_3m

        if len(greater_vals)!=0: 
            mask_conf = greater_vals < threshold
            mask_anti = greater_vals >= threshold

            # convert to copying probabilities
            # greater than majority (adjust according to threshold)
            greater_vals[mask_conf] += (greater_vals[mask_conf] / np.sum(greater_vals)) * D_conf
            greater_vals[mask_anti] += (greater_vals[mask_anti] / np.sum(greater_vals)) * D_anti
            # less than majority
            inv_y = np.zeros_like(less_vals, dtype=float)
            mask = np.abs(less_vals) > 1e-5
            inv_y[mask] = 1 / less_vals[mask]

            if np.sum(inv_y) == 0:
                factors = np.zeros_like(inv_y)
            else:
                factors = inv_y / np.sum(inv_y)
        
            less_vals = less_vals - factors * D_conf
            #print(factors)
            #inv_y = 1/less_vals
            #factors = inv_y / np.sum(inv_y)
            #less_vals = less_vals - factors*D_conf
                
        # construct vector of copying probabilities
        probs = np.empty_like(vec)
        probs[greater_idx] = greater_vals
        probs[equal_idx] = equal_vals
        probs[less_idx] = less_vals

        probs = np.where(probs > 1e-5, probs, 0)
        #print(probs)           
        #print(y_t)

        #normalize
        prob_sum = sum(probs)
        probs = probs/prob_sum
        #print(less_vals, equal_vals, greater_vals, probs)

    return probs
    
def one_generation(m, b, c_t, y_t, q_values, n_population):
    """ simulating n_matings number of matings with a given copying probability"""

    q_proportion = q_values/sum(q_values)
    n_assort = np.round(q_proportion*n_population) # number of females that prefer each male morph (inherent prference)
    n_morphpop = np.round(y_t*n_population) # array of population of each morph
    n_tries = []
    y_og = y_t

    ## memory
    window_size = 10 
    delta_y_history = np.zeros(window_size)

    for mating in range(int(n_matings)):
        counter = 1
        # no copying
        if rnd.random()<(1-c_t):
            
            ### Pair choosing
            chosen_fem = choose(m,n_assort, q_proportion)
            chosen_male = choose(m,n_morphpop,y_t)
        
            if chosen_male!=chosen_fem:
                while(chosen_fem!=chosen_male):
                    counter+=1
                    chosen_fem = choose(m,n_assort, q_proportion)
                    chosen_male = choose(m,n_morphpop,y_t)

        
            # getting number of births
            q = q_values[chosen_male] # geting the male morph's quality
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births # adding births to the male morph population 
        
        else: # copying
            
            switching_probs = copying_probabilities(y_t, b, factor, copying_type = copying_type)
            
            ### Pair choosing
            chosen_fem = choose(m,n_assort,q_proportion)
            # normalizing the probabilities, for copying type 1
            probs = switching_probs
            probs[chosen_fem] = 0
            prob_sum = sum(probs)
            probs[chosen_fem] = np.maximum(1-prob_sum,0)
            
            new_pref = np.random.choice(range(m), p=probs)
            chosen_male = choose(m, n_morphpop, y_t)

            if chosen_male!=new_pref:
                while(new_pref!=chosen_male):
                    counter+=1
                    chosen_fem = choose(m,n_assort,q_proportion)
                    # normalizing the probabilities, for copying type 1
                    probs = switching_probs
                    probs[chosen_fem] = 0
                    prob_sum = sum(probs)
                    probs[chosen_fem] = np.maximum(1-prob_sum,0)
            
                    new_pref = np.random.choice(range(m), p=probs)
                    chosen_male = choose(m,n_morphpop,y_t)

            # getting number of births
            q = q_values[chosen_male] # geting the male morph's quality
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births # adding births to the male morph population 
        
        n_tries.append(counter)
        ### deaths
        n_deaths = n_births # keeping the population constant
        for death in range(n_deaths):
            dying_one = choose(m, n_morphpop, y_t)
            n_morphpop[dying_one] -= 1      
        
        ## update y_t and c_t
        n_population = sum(n_morphpop)
        y_t = n_morphpop/sum(n_morphpop)
        if(1 in y_t):
            break

        # Shift the history array to the left, dropping the oldest value
        diff = np.array(y_t) - np.array(y_og)
        delta_y_history = np.roll(delta_y_history, -1)
        # Place the newest delta_y at the end of the array (now the last position)
        delta_y_history[-1] = diff[-1] # lowest quality morph
        average_delta_y = np.mean(delta_y_history)
    
        ## update copying probability
        if (c_t < 1 and c_t > 0):
            if average_delta_y > 0:
                c_t -= 0.005 
            else:
                c_t += 0.005
        
            #q_avg = np.average(q_values, weights=y_t)
            #q_avg = np.mean(q_values)

            ## 2 & 3
            #if q_values[np.argmax(diff)]>q_avg: dc_t = 0.01 #*np.max(diff)
            #else: dc_t = -0.01 #*np.max(diff)

            ## 1
            #dc_t = np.where(q_values > q_avg, 0.001 , -0.001)
            #dc_t = dc_t[np.argmax(y_t)]

            ## 4
            #if diff[1]>0: dc_t = -0.001
            #else: dc_t = 0.001
            #c_t += dc_t
    #print(np.mean(n_tries))

    return y_t, c_t

def one_run(n_population, q_values, m, T, b, c, y_0):
    c_t = c
    y_t = y_0
    y_hist = [y_0]
    c_hist = [c]
    for t in range(T):
        if(1 not in y_t):
            y_t, c_t = one_generation(m, b, c_t, y_t, q_values, n_population) # updating frequencies after one generation
            y_hist.append(y_t)
            c_hist.append(c_t)
        else:
            y_hist.append(y_t)
            c_hist.append(c_t)

    return y_hist, c_hist

def simulation(n_population, n_runs, m, T, b, c, y_0, q_values, out_dir):
        
    all_runs = []
    all_c = []
    for run in range(n_runs):
        y_hist, c_hist = one_run(n_population, q_values, m, T, b, c, y_0)
        all_runs.append(y_hist)
        all_c.append(c_hist)

    df = pd.DataFrame(all_runs)
    df_c = pd.DataFrame(all_c)
    if m==3:
        df.to_pickle(f'{out_dir}/{float(y_0[0]):.2f}_{float(y_0[2]):.2f}_{float(c):.2f}.pkl')
        df_c.to_pickle(f'{out_dir}/{float(y_0[0]):.2f}_{float(y_0[2]):.2f}_{float(c):.2f}_gamma.pkl')
    if m==2:
        df.to_pickle(f'{out_dir}/{float(y_0[1]):.2f}_{float(c):.2f}.pkl')
        df_c.to_pickle(f'{out_dir}/{float(y_0[1]):.2f}_{float(c):.2f}_gamma.pkl')
    
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