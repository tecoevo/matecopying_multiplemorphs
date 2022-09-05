import math
from matplotlib import pyplot as plt
import ternary
import numpy as np
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import odeint

# plot time heatmap? Plots basins of attraction heatmap if false.
time_heatmap = False 
# resolution of heatmap; the number of points on each side of the simplex
scale = 100

## parameters
# factor of asymmetry
alpha = 0.2
# extent of conformism
beta = 20
# probability of copying
gamma = 0.0
# qualities of male morphs
Q = np.array([2, 2.1, 2.3])
# number of male morphs
m = 3
# fraction of matings in the training period
u = 0.01 

# define colormap to visualize time to fixation
norm = matplotlib.colors.Normalize(vmin=100,vmax=1001)

# define colors to visualize the basins of attraction
blue = '#016fb9'
red = '#e07a5f'
yellow = '#fcec52'
#grey = '#272932'
#white = '#F4F5F5'
colors = [blue, yellow, red]

# define functions
def preference(gamma, y):
    ''' output: array of preference factors for all the morphs; input: gamma, frequency vector'''
    
    q_sum = sum(Q)
    C = 1/(1+ np.exp(-beta*(y-0.5))) # copying function
    P = ((1-gamma)*Q/q_sum) + (gamma*alpha*C)
    return P

def y_dot(y, t):
    ''' output: array of ydots for each male morph population; computes preference factors first using the parameters'''
    
    P = preference(gamma, y) # copying period
    P0 = preference(0, y) # training period

    ydot = y*(u*(P0*Q - np.dot(P0*Q, y)) + (1-u)*(P*Q - np.dot(P*Q, y)))
    return ydot

def solve(x, y, z):
    ''' output: trajectory (array of points) corresponding to the given initial condition [x,y,z]; time is fixed (1000 time steps)'''

    # converting into coordinates
    y_arr = np.array([x,y,z])/100
    # time points
    t = np.arange(0, 1000, 0.1)
    # generating trajectory by integrating
    trajectories = odeint(y_dot, y_arr, t)
    
    return trajectories

def generate_trajectory_data():
    ''' output: dictionary with keys as initial point and values as trajectories; used for plotting with ternary package'''
    
    from ternary.helpers import simplex_iterator
    initpoints = simplex_iterator(scale)

    d = dict()
    print("Solving...")
    for [i, j, k] in initpoints:
        d[(i, j, k)] = solve(i,j,k)
    return d


## begin simulation
data = generate_trajectory_data()

# initialising figure
figure, tax = ternary.figure(scale=scale)
figure.set_size_inches(5, 6)

print("Plotting...")

hm_data = dict()
for key in data.keys():
    traj = data[key]
    if time_heatmap:
        # determining the dominant morph at each time point in the trajectory
        maxes = [max(pt) for pt in traj]
        # time till fixation (dominant morph frequency > 0.99)
        time = min(999, np.argmax(np.array(maxes)>0.99))
        hm_data[key] = time
        if key==sim_point:
            print(time)
    else:
        # dominant morph at the end of the trajectory
        last = np.argmax(np.array(traj[-1]))
        hm_data[key] = last
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))


# plotting he heatmap
if time_heatmap:
    tax.heatmap(hm_data,  style="hexagonal", cmap='jet', use_rgba=False, vmin=100, vmax=1000, cbarlabel="Time to fixation", cb_kwargs=dict(ticks=[100,1000], orientation="horizontal"))
else:
    tax.heatmap(hm_data, style="hexagonal", cmap=cmap, cbarlabel="Basins of Attraction", cticklabels=["Morph 1","Morph 2", "Morph 3"], cb_kwargs=dict(ticks=[1/3,1,5/3], orientation="horizontal"))
   
# adding plot stuff
tax.boundary()
tax.scatter([scale*np.array([1/3,1/3,1/3])], marker='+', c='k', s=10, linewidth=0.5) # center
#tax.scatter([scale*np.array([1/2,1/3,1/6])], c='k', s=20)
#tax.scatter([scale*np.array([0.2,0.2,0.6])], c='k', s=20)
#tax.scatter([scale*np.array([0.4,0.5,0.1])], c='k', s=20)
tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()
#tax.right_corner_label("Morph 1")
#tax.top_corner_label("Morph 2")
#tax.left_corner_label("Morph 3")
#tax.ticks(axis='lbr', ticks=["%.1f" % (2.0 * i / 10) for i in range(6)], offset=0.02)

if time_heatmap:
    tax.savefig("{c}-time.pdf".format(c=str(gamma)))
else:
    tax.savefig("{c}-boa.pdf".format(c=str(gamma)))
