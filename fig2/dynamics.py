import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,100)

b = -20
u = 0.01
a = 2
c = 0.6
q1 = 2
q2 = 2.3

y = x*(1-x)*(((u/(1-u))+c)*(q1**2 - q2**2)/(q1+q2)**2 - (1-c)*a*q2/(q1+q2) + (1-c)*a/(1+np.exp(-b*(x-0.5))))


ax = plt.gca()
ax.set_aspect(1)

plt.plot(x,y,'k-',linewidth=4)
plt.axhline(color='black',linestyle='--',linewidth=2)
plt.box(False)
plt.xticks([])
plt.yticks([])
#plt.plot(x,x,'k--', linewidth=0.5)
plt.savefig("dynamics-anticonformism.svg", format="svg")
