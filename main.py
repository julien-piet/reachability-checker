from src.generics.reachability import Reachability

from src.zonotopes.solver import Solver
s = Solver(0.02, 6)
r = Reachability(s, "examples/example.json")

print(r.automaton)
