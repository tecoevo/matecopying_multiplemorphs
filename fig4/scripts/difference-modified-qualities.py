import math
from matplotlib import pyplot as plt
import ternary
import numpy as np
import matplotlib

# resolution of the heatmap; number of points on one side of the simplex
scale = 10

## parameters
# factor of asymmetry
alpha = 0.2
# extent of conformism
beta = 20
# probability of copying
gamma = 0.5
# qualities of male morphs
Q = np.array([2, 2.1, 2.3])
# number of male morphs
m = 3
# fraction of matings in the training period
u = 0.01 

# defining functions
def preference(gamma, y):
    ''' output: array of preference factors for all the morphs; input: gamma, frequency vector'''
    
    q_sum = sum(Q)
    C = 1/(1+ np.exp(-beta*(y-0.5))) # copying function
    P = ((1-gamma)*Q/q_sum) + (gamma*alpha*C)
    return P

def pref_difference(i,j,k):
    ''' output: difference in the preference-scaled quality values of the 3 morphs (for a given value of gamma)'''

    # converting into coordinates on simplex
    i = i/100
    j = j/100
    k = k/100
    v = np.array([i,j,k])

    pref = preference(gamma, v)
    # converting qualities into vectors 
    vector1 = np.array([0,1])*(Q[0]*pref[0])
    vector2 = np.array([np.sqrt(3)/2, -0.5])*(Q[1]*pref[1])
    vector3 = np.array([-np.sqrt(3)/2, -0.5])*(Q[2]*pref[2])

    return np.linalg.norm(vector1 + vector2 + vector3)
    
def generate_heatmap_data():
    '''output: dictionary with keys as initial point and values as trajectories; used for plotting with ternary package'''

    from ternary.helpers import simplex_iterator
    initpoints = simplex_iterator(scale)

    d = dict()
    print("Solving...")
    for [i, j, k] in initpoints:
        d[(i, j, k)] = pref_difference(i,j,k)
    return d

# generate and plot data
data = generate_heatmap_data()

# initialising figure
figure, tax = ternary.figure(scale=scale)
figure.set_size_inches(5, 6)

print("Plotting...")

# discrete colormap 
cmap = plt.cm.get_cmap('Greys_r', 15)
# plotting the heatmap
tax.heatmap(data,  style="hexagonal", cmap=cmap, use_rgba=False, vmin=0, vmax=0.5,  cbarlabel="Absolute difference in modified qualities", cb_kwargs=dict(ticks=[0.0,0.5], orientation="horizontal"))

# adding plot stuff
tax.boundary()
tax.right_corner_label("Morph 1")
tax.top_corner_label("Morph 2")
tax.left_corner_label("Morph 3")
#tax.ticks(axis='lbr', ticks=["%.1f" % (2.0 * i / 10) for i in range(6)], offset=0.02)
tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()

tax.savefig("{c}-pref.pdf".format(c=str(gamma)))
