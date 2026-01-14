import numpy as np

def generate_grid(n):
    '''Helper function to generate the set of initial points for simulations with 3 morphs'''
    grid = []
    for x in np.linspace(0, 1, n):
        for y in np.linspace(0, 1 - x, n - int(x * n)):
            z = 1 - x - y
            grid.append((x, y, z))
    return grid


### Set parameter values

# Number of male morphs
m = 2 

# beta, the extent of (anti)conformity for Type I copying
# b > 1: conformity
# b = 1: frequency-proportional copying
# b < -1: strong anticonformity
# 0 < b < 1: weak anticonformity
b = 2

# f, the extent of (anti)confrmity for Type II copying
# 0 < f < 1: conformity
# -1 < f < 0: anticonformity
# f = 0: frequency-proportional copying
f = 0.83

# Type of mate copying to determine switching probabilities 
# 1: Type I copying function
# 2: Type II (purely conformist or anticonformist) copying function 
# 3: Type II (mixed conformist and anticonformist) copying function 
copying_type = 3

# Threshold for conformity in Type II mixed copying functions
# For 2 morphs: threshold_2m
# For 3 morphs: threshold_3m
threshold_2m = 0.7
threshold_3m = 0.5

# Strong anticonformity (True/False)
# Only for m=2
strong_anticonformity = False 

# Male population size, total population is twice this number
n_population = 100 
# Number of matings in one generation
n_matings = 100 
# Number of independent runs
n_runs = 50 
# Number of generations
T = 25 

### Set initial conditions

if m==2:
    y_range = np.around(np.arange(0,1.001,0.01), 2) # Initial male morph frequencies
    q_values = np.array([3,2]) # Mate morph qualities
    c_range = np.around(np.arange(0,1.001,0.01), 2) # Mate copying probability (gamma)
if m==3:
    y_range = generate_grid(50) # Initial male morph frequencies
    q_values = np.array([3,2.5,2]) # Mate morph qualities
    c_range = np.around(np.arange(0,1.001,0.05), 2) # Mate copying probability (gamma)
