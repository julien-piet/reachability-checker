class Solver:
    '''
        Interface for an automaton solver.
        Should implement a solve method such as follow :
        Parameters:
            - eq: equation to solve, of type given by get_equation_type()
            - init: set of initial values, of type given by get_set_type()
            - guards: array of possible guards, of type given by get_guard_type()
        Returns:
            - set of initial values for next state
            - index of the guard taken
    '''

    def solve(self, eq, init, guards):
        pass
