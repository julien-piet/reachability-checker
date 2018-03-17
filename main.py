from src.generics.reachability import Reachability
import matplotlib.pyplot as plt

from src.zonotopes.solver import Solver
s = Solver(0.02, 6)
r = Reachability(s, "examples/zono.json")

plt.ion()
plt.draw()
plt.pause(0.001)

r.step(2)
r.step(2)

plt.ioff()
input("Press enter to quit...")
plt.close()
