import math
import numpy as np
import itertools
import pandas as pd

# parameters

alpha = 0.2 # factor of asymmetry
beta = 20 # extent of conformism
Q = np.array([2, 2.3]) # quality vector
m = 2 # number of morphs
u = 0.01 # fraction of matings in the training period

def preference(gamma, y):
    ''' output: array of preference factors for all the morphs; input: gamma, frequency vector'''
    
    q_sum = sum(Q)
    C = 1/(1+ np.exp(-beta*(y-0.5))) # copying function
    P = ((1-gamma)*Q/q_sum) + (gamma*alpha*C)
    return P

def generate_heatmap_data():
    ''' output: dictionary of data (difference in modified qualities) with key values as coordinates on the 2D plane'''
    
    array = np.arange(0,1.01,0.01)
    initpoints = list(itertools.product(array,array)) #coordinates
    
    d = dict()
    for g, i in initpoints: # g=gamma, i=initial y-value
        d[(g,i)] = abs(Q[0]*preference(g, np.array([i,1-i]))[0] - Q[1]*preference(g, np.array([i,1-i]))[1])
        
    return d

data = generate_heatmap_data()
ser = pd.Series(list(data.values()), index=pd.MultiIndex.from_tuples(data.keys()))
df = ser.unstack().fillna(0)

filename = "fig3b_data" 
df.to_pickle(filename+".pkl")
df.to_csv(filename+".csv")
