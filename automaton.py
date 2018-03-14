import numpy as np
from equation import Equation

class Automaton:
    def __init__(self):
        self.nodes = {}
        self.init = ""
        self.x_range = None 

    def set_init_interval(self, inter):
        self.x_range = inter

    def set_init_node(self, node_name):
        self.init = self.nodes[node_name]
        
    def __str__(self):
        rtn = " Automaton : \n"
        for name in self.nodes:
            rtn += str(self.nodes[name])
        rtn += " Init : " + str(self.init)
        return rtn
        

class Link:
    def __init__(self, auto, src, dst):
        self.src = src
        self.src.links.append(self)
        self.dst = dst
        self.guard = None
        self.update = None
        
    def set_update(self, update):
        self.update = update

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
        