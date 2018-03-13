'''
    parser.py - Parse a json to create an automaton object.
'''
import json
import re
import numpy as np

import automaton as am
from interval import Interval
from equation import Equation

def parseFile(filename):
    '''
       Parse a json file to generate an automaton
       For syntax, see example.json
    '''

    data = json.load(open(filename))
    auto = am.Automaton()

    for name in data['nodes']:
        node = am.Node(auto, name)
        if 'eqs' in data['nodes'][name]:
            node.set_equation(parse_system(data['nodes'][name]['eqs']))

    for link in data['links']:
        src = auto.nodes[link['src']]
        dst = auto.nodes[link['dst']]
        link = am.Link(auto, src, dst)

    auto.set_init_interval(parse_init_interval(data['init']))
    auto.set_init_node(data['entry'])

    return auto



eq_name_regex = re.compile("([a-zA-Z]+)' *= *")
eq_params_regex = re.compile("((?:\+|-)? *[0-9]+(?:\.[0-9]*)?) *([a-zA-Z]+)")
eq_bias_regex = re.compile("(?:(±) *([0-9]+(?:\.[0-9]*)?))|(?:(\+|-) *\[ *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *; *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *\])")

def parse_system(eqs):
    '''
        Parse an array of lin dif eq, and create the corresponding n-dimensional equation X' = AX + O(E)
    '''
    system_vars, system = [], {}
    for eq in eqs:
        e = parse_equation(eq)
        if e['var'] not in system_vars:
            system_vars.append(e['var'])
        for v in e['params']:
            if v not in system_vars:
                system_vars.append(v)
        system[e['var']] = e

    system_mat, system_error = [], []
    for v in system_vars:
        eq_line = []
        for p in system_vars:
            if p in system[v]['params']:
                eq_line.append(system[v]['params'][p])
            else:
                eq_line.append(0)
        system_mat.append(eq_line)
        system_error.append(system[v]['error'].toList())

    return Equation(np.transpose(np.array(system_vars)), np.array(system_mat), np.array(system_error))

def parse_equation(eq):
    '''
        Parse a linear differential equation of the form "x = Ax + By + ... (+[a;b] | ±c)"
    '''
    # First, find the variable is equation acts on
    m = re.match(eq_name_regex, eq)
    if not m:
        raise Exception('Invalid equation syntax', eq)
    var_name = m.group(1)

    # Then determine linear parameters
    params = {}
    for mm in re.finditer(eq_params_regex, eq):
        params[mm.group(2)] = float(mm.group(1).replace(' ', ''))

    # Finally, find error
    error = Interval()
    m = re.search(eq_bias_regex, eq)
    if m:
        if m.group(1) != None:
            error.fromMax(float(m.group(2)))
        else:
            error.fromInterval(m.group(3), float(m.group(4)), float(m.group(5)))    
    
    return {'var': var_name, 'params': params, 'error': error}

def parse_update_system(eqs):
    '''
        Parse and array of update, and create the corresponding n-dimensional guard X = AX + B
    '''
    pass

def parse_guard_system(eqs):
    '''
        Parse an array of guard, and create the corresponding n-dimensional guard X >= AX + B
    '''
    pass


init_interval_regex = re.compile("([a-zA-Z]+) *= *\[ *((?:\+|-)? *[[0-9]+(?:\.[0-9]*)?) *; *((?:\+|-)? *[[0-9]+(?:\.[0-9]*)?) *\] *")
def parse_init_interval(eqs):
    intervals = {}
    for eq in eqs:
        m = re.match(init_interval_regex, eq)
        if not m:
            raise Exception('Invalid initial interval syntax', eq)
        intervals[m.group(1)] = Interval()
        intervals[m.group(1)].fromInterval('+', m.group(2), m.group(3))
    return intervals

########################################################################

parseFile('example.json')

    

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
