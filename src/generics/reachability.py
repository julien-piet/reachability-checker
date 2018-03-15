from . import parser

class State:
    def __init__(self, node, values):
        self.node = node
        self.values = values

class Reachability:

    def __init__(self, solver, filename=None):
        self.solver = solver
        self.automaton = None 
        if not filename is None:
            self.load_file(filename)

    def load_file(self, filename):
        eq_type = self.solver.get_equation_type()
        set_type = self.solver.get_set_type()
        guard_type = self.solver.get_guard_type()
        update_type = self.solver.get_update_type()
        self.automaton = parser.parse_file(filename, eq_type, set_type, guard_type, update_type) 

    def step(self):
        if self.states is None:
            assert not self.automaton is None, "Trying to start algorithm on empty automaton"
            self.states = [State(self.automaton.init_node, self.automaton.init_value)]
            return

        next_states = []
        for state in self.states:
            pass
