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
# Extent of (anti)conformity for Type I copying (beta)
b = -2 
# Modulates the extent of conformity/anticonfrmity, higher means less conformity
factor = 3 
# Type of mate copying to determine switching probabilities
copying_type = 3

threshold_2m = 0.7
threshold_3m = 0.5
# Type of anticonformity function (BR type or strong)
BR_type = False 


# Male population size, total population is twice this number
n_population = 100 
# Number of matings in one generation
n_matings = 100 
# Number of independent runs
n_runs = 50 
# Number of generations
T = 50 




if m==2:
    y_range = np.around(np.arange(0,1.001,0.01), 2)
    q_values = np.array([3,2])
    c_range = np.around(np.arange(0,1.001,0.01), 2) # mate copying probability (gamma)
if m==3:
    y_range = generate_grid(50)
    q_values = np.array([3,2.5,2])
    c_range = np.around(np.arange(0,1.001,0.05), 2) # mate copying probability (gamma)
