from sympy import *
import numpy as np


a = 0.2 # factor of asymmetry
b = 20 # extent of conformism
c = 0.2 # probability of copying
u = 0.01 # fraction of matings in training period
q1 = 2 # quality of the lower quality morph
q2 = 2.3 # quality of the higher quality morph

# copying function
def C(x):
    return 1/(1+exp(-b*(x-0.5)))

# frequency of the lower quality morph
y = Symbol('y')

# preference factors
P1 = (1-c)*q1/(q1+q2) + c*a*C(y)
P2 = (1-c)*q2/(q1+q2) + c*a*(1-C(y))

# ydot = f
f = y*(u*((q1**2 - y*(q1**2) - (1-y)*(q2**2))/(q1+q2)) + (1-u)*(P1*q1 - P1*q1*y - P2*q2*(1-y)))

dfdy = diff(f,y)
dfdy = lambdify(y, dfdy)
f = lambdify(y, f)

# check ydot at point z
z = 0
print("ydot at y={z} is {f}".format(z=z, f=f(z)))

# check the stability of fixed point p
p = 0
print("f' at y={p} is {s}".format(p=p, s=dfdy(p)))
