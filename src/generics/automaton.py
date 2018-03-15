# Automaton main class

import numpy as np
from .equation import Equation

class Automaton:

    def __init__(self):
        self.nodes = {}
        self.init_node = None
        self.init_value = None
        self.vars = []
         
        self.solver = None
        self.range = None
        
        self.toExplore = None
        self.explored = None
        
        self.T = None
        self.d = None
        self.N = None

    def set_init_interval(self, inter):
        self.init_value = inter

    def set_variables(self, variables):
        self.variables = variables
        
    def set_init_node(self, node_name):
        self.init = self.nodes[node_name]
        
    def set_solver(self, solv):
        self.solver = solv
        
    def set_limits(self, T, N, d):
        self.T = T
        self.d = d
        self.N = N

    def __str__(self):
        rtn = " Automaton : \n"
        for name in self.nodes:
            rtn += str(self.nodes[name])
        rtn += " Init : " + str(self.init)
        return rtn
        
    # def step():
        # '''
            # Explores current set of nodes
        # '''
        # if self.T is None:
            # print("You need to define serach bounds")
            # return
        # if self.solver is None:
            # print("You need to define a solver")
            # return
            
        # if self.range is None:
            # self.explored = []
            # self.toExplore = [(self.init,self.x_range, 0)]
            # self.range = {}
            
        # current_set = self.toExplore[:]
        # self.toExplore = []
        
        # for val in current_set:
            
            # if val[2] > N:
                # continue
            # # Explore this node
            # solv = self.solver(self.nodes[val[0]], val[1], T, d)
            # solv.run()
            # if not solv.range:
                # continue
                
            # #We have a solution
            # if not val[0] in self.explored:
                # self.explored.append(val[0])
            # if val[0] in self.range:
                # self.range[val[0]] += solv.range
            # else:
                # self.range[val[0]] = solv.range
                
            # #Intersect with guards
            # start = []
            # for set in solv.range:
                # for l in self.nodes[val[0]].links:
                    # new_zone = l.intersect_guard(set)
                    # if not new_zone:
                        # continue
                    # new_zone = l.update(new_zone)
                    # if l.end in starts:
                        # starts[l.end].append(new_zone)
                    # else:
                        # starts[l.end] = [new_zone]
            
            # #Setup exploration for next iteration
            # for nd in starts:
                # for zone in starts[nd]:
                    # self.toExplore.append((nd, zone, val[2]+1))
                    

class Link:
    def __init__(self, auto, src, dst):
        self.src = src
        self.src.links.append(self)
        self.dst = dst
        self.guard = None
        self.update = None
        
    def set_update(self, update):
        self.update = update
        
    def update(self, zone):
        return self.update.update(zone)
        
    def intersect_guard(self, zone):
        return self.guard.intersect(zone)

    def set_guard(self, guard):
        self.guard = guard
        
    def __str__(self):
        rtn = self.src.name + " -> " + self.dst.name
        rtn += "\n Guard : " + str(self.guard)
        rtn += "\n Update : " + str(self.update)
        return rtn
        

class Node:
    def __init__(self, auto, name):
        self.name = name
        self.auto = auto
        self.auto.nodes[name] = self
        self.links = []
        self.equation = None 
        self.guard = None

    def set_equation(self, equation):
        self.equation = equation

    def set_guard(self, guard):
        self.guard = guard
        
    def __str__(self):
        rtn = " Node " + self.name + "\n"
        rtn += str(self.equation) + "\n"
        rtn += "Guard : " + str(self.guard) + "\n"
        for l in self.links:
            rtn += str(l) + "\n"
        return rtn
        
