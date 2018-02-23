'''
    parser.py - Parse a json to create an automaton object. Only for 1D for the moment
    To incorporate time : add another variable, with t_dot = 1
    For ND, should build systems by itself (TODO)
    FIXME : add support for decimal numbers
'''
import json
import re
import numpy as np

def parse(filename="example.json"):
    '''
        Actual parser function
        For syntax see example.json
    '''
    data = json.load(open(filename))
    am = automaton()
    
    for node_name in data["nodes"]:
        am.nodes[node_name] = (node(equation_translator(data["nodes"][node_name]["equation"]),guard_translator(data["nodes"][node_name]["guard"]))
    for link in data["links"]:
        new_link = link(am.nodes[link["src"]], am.nodes[link["dst"]],guard_translator(link["guard"]),update_translator(link["update"]))
        if new_link.src in am.links:
            am.links[new_link.src].append(new_link)
        else:
            am.links[new_link.src] = [new_link]
            
    am.init_node(data["init"]["node"])
    am.init_interval(interval_translator(data["init"]["interval"]))
    
    
def interval_translator(input):
    '''
        Takes an interval [A,B], outputs the interval object
    '''
    regex = re.compile("^\[(?<A>(?<Asign>-|)\d*),(?<B>(?<Bsign>-|)\d*)\]$")
    parts = re.match(regex, input)
    if not parts:
        print("Unrecognized interval")
        return interval(1,0)
    A = int(match.group("A") or 0)
    B = int(match.group("B") or 0)
    Asign = match.group("Asign") or "+"
    Bsign = match.group("Bsign") or "+"
    if Asign == "-":
        A = -A
    if Bsign == "-":
        B = -B

    return interval(A,B)
    
    
def equation_translator(input):
    '''
        Takes an equation under format 'x_dot = Ax + B ± C' and returns an equation object
    '''
    regex = re.compile("^x_dot[ ]*=[ ]*(((?<A>.*?)x)|)[ ]*(?<sign>[\+-]|)[ ]*(?<B>.*?)[ ]*(±(?<C>.*?)|)$")
    parts = re.match(regex, input) 
    if not parts:
        print("Unrecognized equation")
        return equation(0,0,0)
    A = int(match.group("A") or 0)
    B = int(match.group("B") or 0)
    C = int(match.group("C") or 0)
    sign = match.group("sign") or "+"
    
    if sign == "-":
        B = -B
    
    return equation(A, B, C)
    

def guard_translator(input):
    '''
        Takes a guard under format 'Ax + B >=< 0' and returns a guard object
    '''
    regex = re.compile("^((?<A>(?<Asign>-|)\d*)x|)\s*(?<Bsign>[\+-]|)\s*(?<B>\d*|)\s*(?<cmp><|>|<=|>=|=|!=)\s*0$")
    parts = re.match(regex, input) 
    if not parts:
        print("Unrecognized guard")
        return guard(0,0,"!=")
    A = int(match.group("A") or 0)
    B = int(match.group("B") or 0)
    S = match.group("cmp")
    Asign = match.group("Asign") or "+"
    Bsign = match.group("Bsign") or "+"
    
    if Asign == "-":
        A = -A
    if Bsign == "-":
        B = -B
    
    return guard(A,B,S)
    
    
def update_translator(input):
    '''
        Takes an update under format 'X = AX + B' and returns an update object
    '''
    regex = re.compile("^x\s*=\s*((?<A>(?<Asign>-|)\d*)x|)\s*(?<Bsign>[\+-]|)\s*(?<B>\d*|)$")
    parts = re.match(regex, input) 
    if not parts:
        print("Unrecognized update")
        return update(1,0)
    A = int(match.group("A") or 0)
    B = int(match.group("B") or 0)
    Asign = match.group("Asign") or "+"
    Bsign = match.group("Bsign") or "+"
    
    if Asign == "-":
        A = -A
    if Bsign == "-":
        B = -B
    
    return update(A,B)
    
'''
    Usefull classes. FIXME: Will be put in seperate file for final version
'''
    
class automaton:
    '''
        The automaton class is a collection of nodes and links
    '''
    
    def __init__(self):
        self.nodes = {}
        self.links = {}
        self.init = ""
        self.x_range = interval(0,0)
        
    def init_interval(inter):
        self.x_range = inter
    
    def init_node(name):
        self.init = name


class link:
    '''
        Link class
    '''
    
    def __init__(self, src, dst, guard, update):
        self.dst = dst
        self.src = src
        self.guard = guard
        self.update = update
        
    def take(self, X):
        if self.guard.satisfied(X):
            self.update.reset(X)
            return dst
        return src
        

class node:
    '''
        Automaton node
    '''
    
    def __init__(self, equation, guard):
        self.equation = equation
        self.guard = guard


class guard:
    '''
        Node or link guards
        AX + B >=< 0
        
        S in ["<", ">", "=", "!=", ">=", "<="]
    '''
    
    def __init__(self, A, B, S):
        self.A = A
        self.B = B
        self.S = S
        
    def satisfied(self, X):
        value = self.A * X + self.B
        if S == "<":
            return value < 0
        if S == ">":
            return value > 0
        if S == "=":
            return value == 0
        if S == "!=":
            return value != 0
        if S == ">=":
            return value >= 0
        if S == "<=":
            return value <= 0
        print("Unrecognized value of S : " + self.S)
        raise ValueError


class equation:
    '''
        First order differential equation with noize
        X. = AX + B + U(t)
        U = max|U(t)|
    '''
    
    def __init__(self, A, B, U):
        self.B = B
        self.A = A
        self.U = U

class update:
    '''
        Update class
        X = AX + B
    '''
    
    def __init__(self, A, B):
        self.A = A
        self.B = B
        
    
    def reset(self, X):
        return self.A*X + self.B