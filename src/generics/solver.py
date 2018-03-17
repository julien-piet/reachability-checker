from .transition import Transition

class Solver:
    '''
        Interface for an automaton solver.
        Should implement a solve method such as follow :
        Parameters:
            - eq: equation to solve, of type given by get_equation_type()
            - init: set of initial values, of type given by get_set_type()
            - links: array of possible links, with link.guard of type given by get_guard_type()
        Returns:
            - a list of Transitions
    '''

    def solve(self, eq, init, guards):
        pass
