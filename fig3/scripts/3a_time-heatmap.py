import math
import numpy as np
import matplotlib
from matplotlib.colors import LogNorm, Normalize
from scipy.integrate import odeint
import itertools
import pandas as pd

# parameters
alpha = 0.2 # factor of asymmetry
beta = 20 # extent of conformism
Q = np.array([2, 2.3]) # quality vector
m = 2 # number of morphs
u = 0.01 # fraction of matings in the training period

# to normalize time values
norm = matplotlib.colors.Normalize(vmin=0,vmax=1000)

def preference(gamma, y):
    ''' output: array of preference factors for all the morphs; input: gamma, frequency vector'''
    
    q_sum = sum(Q)
    C = 1/(1+ np.exp(-beta*(y-0.5))) # copying function
    P = ((1-gamma)*Q/q_sum) + (gamma*alpha*C)
    return P

def y_dot(y, t, gamma):
    ''' output: array of ydots for each male morph population; computes preference factors first using the parameters'''

    P = preference(gamma, y) # copying period
    P0 = preference(0, y) # training period

    ydot = y*(u*(P0*Q - np.dot(P0*Q, y)) + (1-u)*(P*Q - np.dot(P*Q, y)))
    return ydot

def solve(g, x):
    ''' output: trajectory (array of points) corresponding to the given initial condition [x,y]; number of time steps has been specified'''

    y_arr = np.array([x,1-x])
    gamma = g

    t = np.arange(0, 2000, 0.1)
    trajectories = odeint(y_dot, y_arr, t, args=(gamma,))
    
    return trajectories

def generate_trajectory_data():
    '''output: dictionary with keys as a tuple (gamma, initial y) and values as trajectories '''

    array = np.arange(0,1.01,0.01)
    initpoints = list(itertools.product(array,array))
    
    d = dict()
    print("Solving...")
    for g, i in initpoints:
        d[(g, i)] = solve(g,i)
        
    return d


traj_data = generate_trajectory_data() # trajectory data

time_data = dict()
for key in traj_data.keys():
    traj = traj_data[key]
    last = np.argmax(np.array(traj[-1])) # male morph at fixation
    maxes = [max(pt) for pt in traj] # max frequency at each point
    time = min(1000, np.argmax(np.array(maxes)>0.99))
    time_data[key] = norm(time)


# create dataframe from dict and saving it
ser = pd.Series(list(time_data.values()), index=pd.MultiIndex.from_tuples(time_data.keys()))
df = ser.unstack().fillna(0.99)

filename = 'fig3a_data'
df.to_csv(filename+".csv")
df.to_pickle(filename+".pkl")

