import numpy as np
import random as rnd
import pandas as pd

### Import parameters and initial conditions
import matecopying_params as params

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

### Define functions 

def choose(m, array, probabilities):
    """ Chooses an index from 1 to m with given probabilities 
    and makes sure that array[index] is at least 1 """

    chosen = np.random.choice(range(m), p=probabilities)
    if array[chosen]<1:
        while array[chosen]<1:
            chosen = np.random.choice(range(m), p=probabilities)

    return chosen

def compute_D(greater_vals, f, strong_anticonformity=True):
    ''' Computes the (anti)conformity constant D as defined in Type II copying
    given the male morph frequencies and the extent of (anti)conformity '''

    if not (-1 < f < 1):
        raise ValueError(f"f must be between -1 and 1, got {f}")
    
    y_max = np.max(greater_vals)
    D_max = (1/y_max - 1)*sum(greater_vals)
    D_min = -sum(greater_vals)

    if f>0: 
        D = D_max*f
    elif f<0:
        f_abs = abs(f)
        if strong_anticonformity: 
            D = D_min*f_abs
        else: 
            D = -D_max*f_abs
    else:
        D = 0
    
    return D
    
def switching_probabilities(y_t, b, f, copying_type=1, strong_anticonformity=True):
    """ Computes switching probabilities based on the type of mate copying: 
    1: Type I copying function
    2: Type II (purely conformist or anticonformist) copying function 
    3: Type II (mixed conformist and anticonformist) copying function
    Output: array probabilities of switching to a given morph """

    # Set strong_anticonformity to True for more than 2 morphs
    if len(y_t)>2:
        strong_anticonformity = True

    if copying_type==1:
        if -1 < b <= 0:
            raise ValueError(f"b must be < -1 or > 0, got {b}")
        
        if b >= 1: # Conformity 
            # Note: these probabilities are not normalised here, 
            # but will be normalised later for this case
            probs = (y_t**b)/(y_t**b + (1-y_t)**b)    
        else: # Anticonformity
            if b <= -1: 
                probs = np.where((y_t==0) | (y_t==1), y_t, (y_t**b) / (y_t**b + (1-y_t)**b)) 
            else:
                probs = (y_t**b)/(y_t**b + (1-y_t)**b)
            
            # Normalize the probabilities
            prob_sum = sum(probs)
            probs = probs/prob_sum
            
    elif copying_type==2:
        majority = 1/len(y_t)
        vec = np.array(y_t)
            
        # Partition the frequencies
        greater_idx = np.where(vec > majority)[0]
        equal_idx = np.where(vec == majority)[0]
        less_idx = np.where(vec < majority)[0]

        greater_vals = vec[greater_idx]
        equal_vals = vec[equal_idx]
        less_vals = vec[less_idx]

        # Calculate the switching probabilities
        if len(greater_vals)!=0:
            D = compute_D(greater_vals, f, strong_anticonformity=strong_anticonformity)
            
            greater_vals += (greater_vals/sum(greater_vals))*D
            inv_y = np.zeros_like(less_vals, dtype=float)
            mask = np.abs(less_vals) > 0
            inv_y[mask] = 1/less_vals[mask]

            if np.sum(inv_y) == 0:
                factors = np.zeros_like(inv_y)
            else:
                factors = inv_y/np.sum(inv_y)
                
            less_vals -= factors*D
        
        # Construct vector of switching probabilities
        probs = np.empty_like(vec)
        probs[greater_idx] = greater_vals
        probs[equal_idx] = equal_vals
        probs[less_idx] = less_vals
        probs = np.where(probs < 0, 0, probs)

        # Normalize
        prob_sum = sum(probs)
        probs = probs/prob_sum

    elif copying_type==3:
        majority = 1/len(y_t)
        vec = np.array(y_t)

        if len(y_t)==2: 
            threshold = threshold_2m
        else: 
            threshold = threshold_3m

        # Partition the frequencies 
        greater_idx = np.where(vec > majority)[0]
        equal_idx = np.where(vec == majority)[0]
        less_idx = np.where(vec < majority)[0]

        greater_vals = vec[greater_idx]
        equal_vals = vec[equal_idx]
        less_vals = vec[less_idx]
            
        # Calculate the switching probabilities
        if len(greater_vals)!=0: 
            f_abs = abs(f)
            D_conf = compute_D(greater_vals, f_abs, strong_anticonformity=strong_anticonformity)
            D_anti = compute_D(greater_vals, -f_abs, strong_anticonformity=strong_anticonformity)

            mask_conf = greater_vals < threshold
            mask_anti = greater_vals >= threshold

            if len(greater_vals[mask_anti])!=0:
                greater_vals[mask_conf] += (greater_vals[mask_conf] / np.sum(greater_vals)) * D_conf
                greater_vals[mask_anti] += (greater_vals[mask_anti] / np.sum(greater_vals)) * D_anti

                mask_conf2 = less_vals > 1-threshold
                mask_anti2 = less_vals <= 1-threshold
                
                inv_y = np.zeros_like(less_vals, dtype=float)
                mask = np.abs(less_vals) > 0
                inv_y[mask] = 1/less_vals[mask]

                if np.sum(inv_y) == 0:
                    factors = np.zeros_like(inv_y)
                else:
                    factors = inv_y/np.sum(inv_y)
        
                less_vals[mask_anti2] -= factors*D_anti
                less_vals[mask_conf2] -= factors*D_conf
            else:
                greater_vals += (greater_vals/sum(greater_vals))*D_conf
                inv_y = np.zeros_like(less_vals, dtype=float)
                mask = np.abs(less_vals) > 0
                inv_y[mask] = 1/less_vals[mask]

                if np.sum(inv_y) == 0:
                    factors = np.zeros_like(inv_y)
                else:
                    factors = inv_y/np.sum(inv_y)
        
                less_vals = less_vals - factors*D_conf
    
            
        # Construct the vector of switching probabilities
        probs = np.empty_like(vec)
        probs[greater_idx] = greater_vals
        probs[equal_idx] = equal_vals
        probs[less_idx] = less_vals
        probs = np.where(probs > 0, probs, 0)

        # Normalize
        prob_sum = sum(probs)
        probs = probs/prob_sum

    return probs
    
def one_generation(c, y_t, n_matings):
    """ Simulates n_matings number of matings with a given 
    mate copying probability (c), starting with the 
    initial morph frequencies y_t """

    # Compute the number of females that prefer each male morph (inherent prference)
    q_proportion = q_values/sum(q_values)
    n_assort = np.round(q_proportion*n_population) 

    # Population size of each male morph
    n_morphpop = np.round(y_t*n_population) 
    
    # Counts for the number of tries before successful mating
    n_tries = []

    for mating in range(n_matings):
        counter = 1

        # No mate copying with probability (1-c)
        if rnd.random()<(1-c):
            
            # Choose a pair for mating
            chosen_fem = choose(m,n_assort, q_proportion)
            chosen_male = choose(m,n_morphpop,y_t)
        
            if chosen_male!=chosen_fem:
                while(chosen_fem!=chosen_male):
                    counter+=1
                    chosen_fem = choose(m,n_assort, q_proportion)
                    chosen_male = choose(m,n_morphpop,y_t)

            # Add births to the male morph population 
            q = q_values[chosen_male]
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births 
        
        # Mate copying with probability c
        else: 
            switching_probs = switching_probabilities(y_t, b = b, f = f, copying_type = copying_type, strong_anticonformity=strong_anticonformity)
            
            # Choose a pair for mating
            chosen_fem = choose(m,n_assort,q_proportion)
            
            # Normalize the probabilities for Type I conformist copying
            # They depend on the initial peference of the chossing female
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
                    
                    # Normalize the probabilitie for Type I conformist copying
                    probs = switching_probs
                    probs[chosen_fem] = 0
                    prob_sum = sum(probs)
                    probs[chosen_fem] = np.maximum(1-prob_sum,0)
            
                    new_pref = np.random.choice(range(m), p=probs)
                    chosen_male = choose(m,n_morphpop,y_t)

            # Add births to the male morph population
            q = q_values[chosen_male] 
            n_births = np.random.poisson(q,1)[0]
            n_morphpop[chosen_male] += n_births 
        
        n_tries.append(counter)

        ### Deaths, keeping the population at a constant size
        n_deaths = n_births 
        for death in range(n_deaths):
            dying_one = choose(m, n_morphpop, y_t)
            n_morphpop[dying_one] -= 1      
        
        y_t = n_morphpop/sum(n_morphpop)
        
        if(1 in y_t):
            break

    return y_t

def one_run(c, y_0, T):
    """ Simulates T number of generations with a given 
    mate copying probability (c), starting with the 
    initial morph frequencies y_0 """

    y_t = y_0
    y_hist = [y_0]
    for t in range(T):
        if(1 not in y_t):
            y_t = one_generation(c, y_t, n_matings=n_matings) # Update frequencies after one generation
            y_hist.append(y_t)
        else:
            y_hist.append(y_t)

    return y_hist

def simulation(c, y_0, n_runs, out_dir):
    """ Simulates n_runs number of runs with a given 
    mate copying probability (c), starting with the 
    initial morph frequencies y_0. Stores the output
    in out_dir for 2 and 3 morphs """

    all_runs = []
    for run in range(n_runs):
        y_hist = one_run(c, y_0, T=T)
        all_runs.append(y_hist)

    df = pd.DataFrame(all_runs)
    if m==3:
        df.to_pickle(f'{out_dir}/{float(y_0[0]):.2f}_{float(y_0[2]):.2f}_{float(c):.2f}.pkl')
    if m==2:
        df.to_pickle(f'{out_dir}/{float(y_0[1]):.2f}_{float(c):.2f}.pkl')
    
def run_simulations_2m(y,c, out_dir):
    simulation(c, np.array([1-y,y]), n_runs, out_dir=out_dir)
    return f"Output written for y={y} and c={c}"

def run_simulations_3m(y,c, out_dir):
    simulation(c, np.array(y), n_runs, out_dir)
    return f"Output written for y={y} and c={c}"

def fix_fractions(x,y,c,input_dir):
    """ Computes the average outcome of several independent runs for a given
    initial condition """

    f = f'{float(x):.2f}_{float(y):.2f}_{float(c):.2f}.pkl'
    df = pd.read_pickle(f'{input_dir}/{f}')
    last_gen2 = np.sum([df[df.shape[1]-1][i][2] for i in range(df.shape[0])])/50
    last_gen0 = np.sum([df[df.shape[1]-1][i][0] for i in range(df.shape[0])])/50
    last_gen1 = np.sum([df[df.shape[1]-1][i][1] for i in range(df.shape[0])])/50
    return last_gen0, last_gen1, last_gen2