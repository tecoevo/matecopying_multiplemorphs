import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import itertools
import re

allfiles = os.listdir()
datafiles = []

# get the data files from all files in the directory
for f in allfiles:
    match = re.search('simulation_data_.*pkl',f)
    if match!=None:
        datafiles.append(f)

files = sorted(datafiles) 

# range of gamma values
param_values_c = np.around(np.arange(0, 1.01, 0.01), 2)
# lattice of y and gamma values
params_lattice = list(itertools.product(param_values_c, param_values_c))

for f in files:
    coordinate = params_lattice[files.index(f)]
    df = pd.read_pickle(f)

    # average the outcome over all the runs
    last_gen = sum(np.array(df[df.shape[1]-1]))
    ind = np.argmax(last_gen)

    # colour according to outcome
    if ind==0: 
        color = '#588EF2'
    else: 
        color = '#E9B360'

    with plt.style.context('bmh'):
        plt.scatter(coordinate[1], coordinate[0], c=color, marker='s', s=50)

with plt.style.context('bmh'):
    plt.ylabel('Initial frequency of the lower quality morph')
    plt.xlabel('$\gamma$')
    plt.title('Bifurcation diagram (simulation)')

    ticks = np.arange(0,1.1,0.1)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.grid(alpha=0.5)
    
    
u = 0.1 # fraction of matings in the training period
q1 = 2 # lower quality morph
q2 = 2.3 # higher quality morph
a = 0.4 # factor of asymmetry
b = 4 # extent of conformism

# range of gamma values for plotting model output
c_range = np.linspace(0,1,1000)

polymorph = []
bif_pt = 0
is_bif_pt = False

for c in c_range:
    k = -1 + (q1+q2)/(q2 + ((q2-q1)*(u/(1-u) + 1-c)/(a*c)))
    y =  0.5 - (np.log(k)/b) # internal fixed point

    # characterising bifurcation point
    if 1-y>0.0001:
        if is_bif_pt==False:
            print('Bifurcation point: ', c)
            bif_pt = c
            is_bif_pt = True
        
    polymorph.append(y)

after_bif = [x for x in c_range if x>bif_pt]
before_bif = [x for x in c_range if x<=bif_pt]

for i in range(len(before_bif)):
    polymorph[i] = 1

ax = plt.gca()
ax.set_aspect(1)
plt.plot(c_range, polymorph, 'k', alpha=0.5)
plt.xlim([0,1])
plt.ylim([0,1])
plt.savefig('fig-sim.svg', format="svg")

