import numpy as np

def generate_grid(n):
    grid = []
    for x in np.linspace(0, 1, n):
        for y in np.linspace(0, 1 - x, n - int(x * n)):
            z = 1 - x - y
            grid.append((x, y, z))
    return grid

### parameter values
n_population = 100 # male population, total population is twice this number
n_matings = 100 # number of matings in one generation
n_runs = 50 # number of independent runs
m = 2 # number of male morphs
T = 25 # number of generations
b = 2 # extent of conformism (beta)
copying_type = 1 # type of mate copying to determine copying probabilities
factor = 2 # modulates the extent of conformity/anticonfrmity, higher means less conformity
threshold_2m = 0.7
threshold_3m = 0.5

if m==2:
    y_range = np.around(np.arange(0,1.001,0.01), 2)
    q_values = np.array([3,2])
    c_range = np.around(np.arange(0,1.001,0.05), 2) # mate copying probability (gamma)
if m==3:
    y_range = generate_grid(50)
    q_values = np.array([3,2.5,2])
    c_range = np.around(np.arange(0,1.001,0.05), 2) # mate copying probability (gamma)
