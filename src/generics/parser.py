'''
    parser.py - Parse a json to create an automaton object.
'''
import json
import re
import numpy as np

from . import automaton as am
from .interval import Interval


def parse_file(filename, eq_type, set_type, guard_type, update_type):
    '''
       Parse a json file to generate an automaton
       For syntax, see the README 
    '''
    text = open(filename)
    return parse(text, eq_type, set_type, guard_type, update_type)


def parse(text, eq_type, set_type, guard_type, update_type):
    '''
        Same as parse_file, but directly takes file content
    '''

    data = json.load(text)
    auto = am.Automaton()

    auto.set_variables(data['vars'])
    auto.set_init_interval(parse_init_interval(data['init'], data['vars'], set_type))

    for name in data['nodes']:
        node = am.Node(auto, name)
        if 'eqs' in data['nodes'][name]:
            node.set_equation(parse_system(data['nodes'][name]['eqs'], data['vars'], eq_type))
        if 'guard' in data['nodes'][name]:
            node.set_guard(parse_guard_system(data['nodes'][name]['guard'], data['vars'], guard_type))

    for link in data['links']:
        src = auto.nodes[link['src']]
        dst = auto.nodes[link['dst']]
        link_obj = am.Link(auto, src, dst)
        if 'update' in link:
            link_obj.set_update(parse_update_system(link['update'], data['vars'], update_type))
        if 'guard' in link:
            link_obj.set_guard(parse_guard_system(link['guard'], data['vars'], guard_type))
        

    auto.set_init_node(data['entry'])

    return auto

def get_var(variables, v):
    try:
        i = variables.index(v)
        return i
    except ValueError:
        raise Exception("Undeclared variable in equation: " + v)

###############################################################################

eq_name_regex = re.compile("([a-zA-Z]+)' *= *")
eq_params_regex = re.compile("((?:\+|-)? *[0-9]+(?:\.[0-9]*)?)? *([a-zA-Z]+)")
eq_bias_regex = re.compile("(?:(±) *([0-9]+(?:\.[0-9]*)?))|(?:(\+|-) *\[ *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *; *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *\])")

def parse_system(eqs, variables, eq_type):
    '''
        Parse an array of lin dif eq, and create the corresponding n-dimensional equation X' = AX + O(E)
    '''
    A = np.zeros((len(variables), len(variables)))
    B = [[0, 0] for i in range(len(variables))]
    for eq in eqs:
        e = parse_equation(eq)
        i = get_var(variables, e['var'])
        for p in e['params']:
            j = get_var(variables, p)
            A[i][j] = e['params'][p]
        B[i] = e['error'].toList()

    return eq_type.fromSystem(A, np.array(B))

def parse_equation(eq):
    '''
        Parse a linear differential equation of the form "x' = Ax + By + ... (+[a;b] | ±c)"
    '''
    # First, find the variable is equation acts on
    m = re.match(eq_name_regex, eq)
    if not m:
        raise Exception('Invalid equation syntax', eq)
    var_name = m.group(1)

    # Then determine linear parameters
    params = {}
    for mm in re.finditer(eq_params_regex, eq):
        i = float(mm.group(1).replace(' ', '')) if not mm.group(1) is None else 1.0
        params[mm.group(2)] = i

    # Finally, find error
    error = Interval()
    m = re.search(eq_bias_regex, eq)
    if m:
        if m.group(1) != None:
            error.fromMax(float(m.group(2)))
        else:
            error.fromInterval(m.group(3), float(m.group(4)), float(m.group(5)))    

    return {'var': var_name, 'params': params, 'error': error}

###############################################################################

update_name_regex = re.compile("([a-zA-Z]+) *= *")
update_offset_regex = re.compile("([+-]? *[0-9]+(?:\.[0-9]*)?) *$")

def parse_update_system(eqs, variables, update_type):
    '''
        Parse and array of update, and create the corresponding n-dimensional guard X = AX + B
    '''
    A = np.zeros((len(variables), len(variables)))
    B = np.array([0 for i in range(len(variables))])
    for eq in eqs:
        u = parse_update(eq)
        i = get_var(variables, u['var'])
        for p in u['params']:
            j = get_var(variables, p)
            A[i][j] = u['params'][p]
        B[i] = u['offset']

    return update_type.fromSystem(A, B)
    

def parse_update(eq):
    '''
        Parse a single update equation, of the form x = ax + by + ... + c
    '''
    # First, find the variable is equation acts on
    m = re.match(update_name_regex, eq)
    if not m:
        raise Exception('Invalid update syntax', eq)
    var_name = m.group(1)

    # Then determine linear parameters
    params = {}
    for mm in re.finditer(eq_params_regex, eq):
        i = float(mm.group(1).replace(' ', '')) if not mm.group(1) is None else 1.0
        params[mm.group(2)] = i

    # Finally, find offset
    offset = 0
    m = re.search(update_offset_regex, eq)
    if m:
        offset = float(m.group(1).replace(' ',''))  

    return {'var': var_name, 'params': params, 'offset': offset}

###############################################################################

def parse_guard_system(eqs, variables, guard_type):
    '''
        Parse and array of guard, and create the corresponding n-dimensional guard AX + B cmp 0
    '''
    A, B, C = [], [], []
    for eq in eqs:
        g = parse_guard(eq)
        AA = np.zeros(len(variables))
        for p in g['params']:
            j = get_var(variables, p)
            AA[j] = g['params'][p]
        A.append(AA)
        B.append(g['offset'])
        C.append(g['comp'])

    return guard_type.fromSystem(np.array(A), np.array(B), np.array(C))
    
comparator_regex = re.compile("(<|>|=)")
guard_offset_regex = re.compile("([+-] *[0-9]+(?:\.[0-9]*)?) *(<|>|=)")

def parse_guard(eq):
    '''
        Parse a single guard equation, of the form ax + by + ... + c cmp 0
    '''
    # First, find the comparator 
    m = re.search(comparator_regex, eq)
    if not m:
        raise Exception('Invalid guard syntax', eq)
    comp = m.group(1)

    # Then determine linear parameters
    params = {}
    for mm in re.finditer(eq_params_regex, eq):
        i = float(mm.group(1).replace(' ', '')) if not mm.group(1) is None else 1.0
        params[mm.group(2)] = i

    # Finally, find offset
    offset = 0
    m = re.search(guard_offset_regex, eq)
    if m:
        offset = float(m.group(1).replace(' ',''))  
    
    return {'comp': comp, 'params': params, 'offset': offset }

###############################################################################

init_interval_regex = re.compile("([a-zA-Z]+) *= *\[ *((?:\+|-)? *[[0-9]+(?:\.[0-9]*)?) *; *((?:\+|-)? *[[0-9]+(?:\.[0-9]*)?) *\] *")

def parse_init_interval(eqs, variables, set_type):
    I = np.zeros((len(variables), 2))
    for eq in eqs:
        m = re.match(init_interval_regex, eq)
        if not m:
            raise Exception('Invalid initial interval syntax', eq)
        i = get_var(variables, m.group(1))
        intv = Interval()
        intv.fromInterval('+', m.group(2), m.group(3))
        I[i][0] = intv.a
        I[i][1] = intv.b

    return set_type.fromIntervals(I)

