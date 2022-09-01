import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,100)

b = -20

# logistic copying function
y = 1/(1+np.exp(-b*(x-0.5))) 


ax = plt.gca()
ax.set_aspect(1)

plt.plot(x,y,'k-',linewidth=4)
plt.plot(x, x, 'k--', linewidth=2)

plt.xticks([])
plt.yticks([])
#plt.plot(x,x,'k--', linewidth=0.5)
plt.savefig("copying-fn-anticonformism.svg", format="svg")
