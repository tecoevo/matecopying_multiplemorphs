import numpy as np
import matplotlib.pyplot as plt

# extent of conformism (beta)
b = 20

# logistic copying function
x = np.linspace(0,1,100)
y = 1/(1+np.exp(-b*(x-0.5)))

fig, ax = plt.subplots()
with plt.style.context('bmh'):
    if b<0:
        plt.plot(x,y,'k--',label=r'$\beta = $'+str(b), linewidth=0.5, alpha=0.6)
    elif b==0:
        plt.plot(x,y, linewidth=5)
    else:
        plt.plot(x,y,linewidth=4)

plt.grid()
plt.xticks([])
plt.yticks([])
fig.savefig("copying-fn.svg", format="svg")
