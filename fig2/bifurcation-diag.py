import numpy as np
import matplotlib.pyplot as plt

# plot what?
conformism = False # anticonformism is false

# parameters
u = 0.01 # fraction of matings in the training period
q1 = 2 # quality of male morph 1
q2 = 2.3 # quality of male morph 2
a = 0.18 # factor of assymetry
if conformism: b = 20 # extent of conformism
else: b = -20

c_arr = np.linspace(0,1,1000) # range of gamma values
Y = [] # values of the internal fixed point for different gamma values
Z = [] # to plot the values of the first fixed point

for c in c_arr:
    k = -1 + (q1+q2)/(q2 + ((q2-q1)*(u/(1-u) + 1-c)/(a*c)))
    y = min(0.5 - (np.log(k)/b), 1)
    if y<1: z = 1
    else: z = 0

    Z.append(z)
    Y.append(y) 

Z = np.array(Z)
c_arr_afterbif = c_arr[np.nonzero(Z)]

Y = np.array(Y)
if conformism: Y[np.isnan(Y)]=1
else: Y[np.isnan(Y)]=0

ax = plt.gca()
ax.set_aspect(1) # set aspect ratio

if conformism:
    plt.plot(c_arr, Y, color='k',linewidth=3, linestyle='--', dashes=(7,3))
    plt.plot(c_arr, [0]*len(c_arr), 'k-', linewidth=3)
    plt.plot(c_arr_afterbif, Z[np.nonzero(Z)], 'k-',linewidth=3)
else:
    plt.plot(c_arr, Y, 'k-',linewidth=3)
    plt.plot(c_arr, [1]*len(c_arr),color='k',linewidth=3, linestyle='--', dashes=(7,3))
    plt.plot(c_arr, [0]*len(c_arr),color='k',linewidth=3, linestyle='--', dashes=(7,3))

plt.xlim([0,1])
plt.ylim([-0.008,1.008])
if conformism: plt.fill_between(c_arr,Y,np.max(Y), color='#E9B360')
else: plt.fill_between(c_arr, [1]*len(c_arr), color='#E9B360')
plt.fill_between(c_arr,Y,[0]*len(c_arr), color='#588EF2')

if conformism: filename="conformism.pdf"
else: filename="anticonformism.pdf"
plt.savefig(filename, format='pdf')
#plt.show()
