from . import parser
import matplotlib.pyplot as plt
import numpy
import random

class State:
    def __init__(self, node, values):
        self.node = node
        self.values = values

    def __str__(self):
        return "To node :" + self.node.name + "; values :" + str(self.values)

class Reachability:

    def __init__(self, solver, filename=None):
        self.solver = solver
        self.automaton = None 
        self.colors = {}
        self.T = 1
        self.states = None
        if not filename is None:
            self.load_file(filename)

    def load_file(self, filename):
        eq_type = self.solver.get_equation_type()
        set_type = self.solver.get_set_type()
        guard_type = self.solver.get_guard_type()
        update_type = self.solver.get_update_type()
        self.automaton = parser.parse_file(filename, eq_type, set_type, guard_type, update_type) 

    def step(self, T):
        if self.states is None:
            assert not self.automaton is None, "Trying to start algorithm on empty automaton"
            self.states = [State(self.automaton.init_node, self.automaton.init_value)]

        next_states = []
        for state in self.states:
            if not state.node.name in self.colors:
                r = random.random()
                b = random.random()
                g = random.random()
                self.colors[state.node.name] = (r, g, b)

            eq = state.node.equation
            val = state.values
            node = state.node
            trs = self.solver.solve(eq, val, node, T, self.colors[state.node.name])
            for tr in trs:
                next_states.append(State(tr.link.dst, tr.values))

        self.states = next_states
