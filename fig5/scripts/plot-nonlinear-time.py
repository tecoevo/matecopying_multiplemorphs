from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
import sys

# define colors to visualize the basins of attraction
blue = '#016fb9'
red = '#e07a5f'
yellow = '#fcec52'
colors = [blue, yellow, red]

file1name = sys.argv[1]
data = pd.read_csv(file1name)
y = [data.index[data['basin']==0], data.index[data['basin']==1], data.index[data['basin']==2]]


switch = False
basin = 2
switchvalue=1
for i in range(999):
    print(i)
    if switch==True:
        break
    else:
        if i in y[0]:
            switch=True
            switchvalue = data['gamma'][i]
            basin = 0
        elif i in y[1]:
            switch=True
            switchvalue = data['gamma'][i]
            basin = 1
print("Plotting...")
with plt.style.context('ggplot'):
    plt.subplot(3,1,1)
    plt.plot(np.array(data['gamma']), np.array(data['time']), 'k', linewidth=2)
    plt.xlim([0,1])
    plt.axvline(switchvalue, color='k', alpha=0.1)
    plt.axvspan(0, switchvalue, facecolor=red)
    plt.axvspan(switchvalue,1, facecolor=colors[basin])
    plt.yscale('log')
    plt.yticks([100, 1000], color="k")
    plt.xticks(color="k")
    plt.minorticks_off()
    plt.grid()
    plt.savefig("pt3.svg", format="svg")
    

## file 2
# file2name = sys.argv[2]
# data = pd.read_csv(file2name)    
# y = [data.index[data['basin']==0], data.index[data['basin']==1], data.index[data['basin']==2]]


# switch = False
# basin = 2
# switchvalue=1
# for i in range(999):
#     print(i)
#     if switch==True:
#         break
#     else:
#         if i in y[0]:
#             switch=True
#             switchvalue = data['gamma'][i]
#             basin = 0
#         elif i in y[1]:
#             switch=True
#             switchvalue = data['gamma'][i]
#             basin = 1
# print("Plotting again...")
# with plt.style.context('ggplot'):
#     plt.subplot(3,1,2)
#     plt.plot(np.array(data['gamma']), np.array(data['time']), 'k', linewidth=2)
#     plt.xlim([0,1])
#     plt.axvline(switchvalue, color='k', alpha=0.1)
#     plt.axvspan(0, switchvalue, facecolor=red)
#     plt.axvspan(switchvalue,1, facecolor=colors[basin])
#     plt.yticks([100, 1500])
#     plt.grid()
#     plt.xticks([])
#     plt.yscale('log')
#     plt.yticks([100, 1000])
#     plt.minorticks_off()
    

# data = pd.read_csv('5.csv')
# y = [data.index[data['basin']==0], data.index[data['basin']==1], data.index[data['basin']==2]]


# switch = False
# basin = 2
# switchvalue=1
# for i in range(999):
#     if switch==True:
#         break
#     else:
#         if i in y[0]:
#             switch=True
#             switchvalue = data['gamma'][i]
#             basin = 0
#         elif i in y[1]:
#             switch=True
#             switchvalue = data['gamma'][i]
#             basin = 1

# with plt.style.context('ggplot'):
#     plt.subplot(3,1,3)
#     plt.plot(np.array(data['gamma']), np.array(data['time']), 'k', linewidth=2)
#     plt.xlim([0,1])
#     plt.axvline(switchvalue, color='k', alpha=0.1)
#     plt.axvspan(0, switchvalue, facecolor=red)
#     plt.axvspan(switchvalue,1, facecolor=colors[basin])
#     plt.yticks([0, 5000])
#     plt.grid()
#     plt.xticks([])
#     plt.yscale('log')
#     plt.yticks([100, 1000])
#     plt.minorticks_off()
    



    
