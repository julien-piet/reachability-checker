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

class Link:
    def __init__(self, auto, src, dst):
        self.src = src
        self.src.links.append(self)
        self.dst = dst

class Node:
    def __init__(self, auto, name):
        self.name = name
        self.auto = auto
        self.auto.nodes[name] = self
        self.links = []
        self.equation = None 

    def set_equation(self, equation):
        self.equation = equation
