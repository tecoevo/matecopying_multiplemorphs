import numpy as np
import matplotlib.pyplot as plt

# extent of conformism, beta
b = -20

# logistic copying function
x = np.linspace(0,1,100)
y = 1/(1+np.exp(-b*(x-0.5))) 

ax = plt.gca()
ax.set_aspect(1)
plt.plot(x,y,'k-',linewidth=4)
plt.plot(x, x, 'k--', linewidth=2)
plt.xticks([])
plt.yticks([])

if b<0:
    plt.savefig("copying-fn-anticonformism.svg", format="svg")
else:
    plt.savefig("copying-fn-conformism.svg", format="svg")
