import math
import numpy as np
from scipy.integrate import odeint
import csv


# point to plot
sim_point = (0.3, 0.6, 0.1) # lowest to hisghest quality

## parameters
alpha = 0.2 # factor of asymmetry 
beta = 20 # extent of conformism
Q = np.array([2, 2.1, 2.3]) # quality vector
m = 3 # number of morphs
u = 0.01 # fraction of matings in the training period

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

def solve(x, y, z, gamma):
    ''' output: trajectory (array of points) corresponding to the given initial condition [x,y,z]; time has been defined (1000 time steps)'''

    y_arr = np.array([x,y,z])

    t = np.arange(0, 1000, 0.1)
    trajectories = odeint(y_dot, y_arr, t, args=(gamma,))
    
    return trajectories

print("Calculating Time...")
data = []

for gamma in np.linspace(0,1,1000):
    traj = solve(sim_point[0], sim_point[1], sim_point[2], gamma)
    maxes = [max(pt) for pt in traj]
    time = np.argmax(np.array(maxes)>0.99)
    last = np.argmax(np.array(traj[-1]))

    data.append({'gamma':gamma, 'time':time, 'basin':last})

print('Writing...')
fields = ['gamma', 'time', 'basin']
with open('l{l}-m{m}-h{h}.csv'.format(l=str(sim_point[0]), m=str(sim_point[1]), h=str(sim_point[2])), 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = fields)
    writer.writeheader()
    writer.writerows(data)



