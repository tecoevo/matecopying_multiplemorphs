import numpy as np
import matplotlib.pyplot as plt

# extent of conformism
b = 20

x = np.linspace(0,1,100)
y1 = 1/(1+np.exp(-b*(x-0.5))) # copying function
y2 = x**2 + (y1*(1-x**2 - (1-x)**2)) # effective conformism function

plt.plot(x, y1, "k--", label="$C(y)$")
plt.plot(x, y2, "k", label="$C_{eff}(y)$ for $s=2$")
plt.xlabel("Frequecy of the male morph ($y$)")
plt.ylabel("Copying function")
plt.legend()
plt.savefig("effective-conformism.pdf", format="pdf")
