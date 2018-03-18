from src.generics.reachability import Reachability
import matplotlib.pyplot as plt
import numpy as np

from src.polygons.continuous_range_calculator import solver
s = solver(0.1)
r = Reachability(s, "examples/example-1D.json")

plt.ion()
plt.draw()
plt.pause(0.001)

r.step(2)
r.step(2)

t = np.linspace(0,2,num=50)
exp = [1.1 * np.exp(i) - 0.1 for i in t]
plt.plot(t,exp,c="b")

plt.ioff()
input("Press enter to quit...")
plt.close()
