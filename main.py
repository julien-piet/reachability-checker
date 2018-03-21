from src.generics.reachability import Reachability
import matplotlib.pyplot as plt
import numpy as np

from src.polygons.continuous_range_calculator import solver
s = solver(0.2)
r = Reachability(s, "examples/thermostat.json")

plt.ion()
plt.draw()
plt.pause(0.001)

r.step(20)
r.step(20)
r.step(20)
r.step(20)

# def exp(alpha, beta, t0, x0, t):
#     return [(beta / alpha) * ( np.exp( alpha * (i - t0) ) - 1 ) + x0 * np.exp( alpha * (i - t0) ) for i in t]
# 
# t = np.linspace(0,5,num=5000)
# plt.plot(t,exp(1,0,0,1,t),c="b")
# plt.plot(t,5*t+1,c="r")
# plt.plot(t,5*t + 5*(1 - np.log(5)), c="g")

plt.ioff()
input("Press enter to quit...")
plt.close()
