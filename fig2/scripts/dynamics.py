import numpy as np
import matplotlib.pyplot as plt

# parameters
u = 0.01 # training period
a = 2 # alpha, factor of asymmetry
c = 0.6 # gamma, probability of copying
q1 = 2 # quality of morph 1
q2 = 2.3 # quality of morph 2

# extent of conformism (conformism if positive, anticonformism if negative)
b = -20 

# dynamics of the male morph population
x = np.linspace(0,1,100)
y = x*(1-x)*(((u/(1-u))+c)*(q1**2 - q2**2)/(q1+q2)**2 - (1-c)*a*q2/(q1+q2) + (1-c)*a/(1+np.exp(-b*(x-0.5))))

# plotting
ax = plt.gca()
ax.set_aspect(1)
plt.plot(x,y,'k-',linewidth=4)
plt.axhline(color='black',linestyle='--',linewidth=2)
plt.box(False)
plt.xticks([])
plt.yticks([])

if b<0:
    plt.savefig("dynamics-anticonformism.svg", format="svg")
else:
    plt.savefig("dynamics-conformism.svg", format="svg")
