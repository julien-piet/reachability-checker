'''
    parser.py - Parse a json to create an automaton object.
'''
import json
import re
import numpy as np

import automaton as am
from interval import Interval
from equation import Equation
from update import Update
from guard import Guard

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
        if 'guard' in data['nodes'][name]:
            node.set_guard(parse_guard_system(data['nodes'][name]['guard']))

    for link in data['links']:
        src = auto.nodes[link['src']]
        dst = auto.nodes[link['dst']]
        link_obj = am.Link(auto, src, dst)
        if 'update' in link:
            link_obj.set_update(parse_update_system(link['update']))
        if 'guard' in data['nodes'][name]:
            link_obj.set_guard(parse_guard_system(link['guard']))
        

    auto.set_init_interval(parse_init_interval(data['init']))
    auto.set_init_node(data['entry'])

    return auto


eq_name_regex = re.compile("([a-zA-Z]+)' *= *")
eq_params_regex = re.compile("((?:\+|-)? *[0-9]+(?:\.[0-9]*)?) *([a-zA-Z]+)")
eq_bias_regex = re.compile("(?:(±) *([0-9]+(?:\.[0-9]*)?))|(?:(\+|-) *\[ *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *; *((?:\+|-)?[0-9]+(?:\.[0-9]*)?) *\])")

update_name_regex = re.compile("([a-zA-Z]+) *= *")
update_offset_regex = re.compile("([+-] *[0-9]+) *$")

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
    system_vars, system = [], {}
    for eq in eqs:
        e = parse_update(eq)
        if e['var'] not in system_vars:
            system_vars.append(e['var'])
        for v in e['params']:
            if v not in system_vars:
                system_vars.append(v)
        system[e['var']] = e

    system_mat, system_offset = [], []
    for v in system_vars:
        eq_line = []
        for p in system_vars:
            if p in system[v]['params']:
                eq_line.append(system[v]['params'][p])
            else:
                eq_line.append(0)
        system_mat.append(eq_line)
        system_offset.append(system[v]['offset'])

    return Update(np.transpose(np.array(system_vars)), np.array(system_mat), np.array(system_offset))
    

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
        params[mm.group(2)] = float(mm.group(1).replace(' ', ''))

    # Finally, find offset
    offset = 0
    m = re.search(update_offset_regex, eq)
    if m:
        offset = float(m.group(1).replace(' ',''))  
    
    return {'var': var_name, 'params': params, 'offset': offset}


def parse_guard_system(eqs):
    '''
        Parse and array of guard, and create the corresponding n-dimensional guard AX + B cmp 0
    '''
    system_vars, system = [], []
    for eq in eqs:
        e = parse_guard(eq)
        for v in e['params']:
            if v not in system_vars:
                system_vars.append(v)
        system.append(e)

    system_mat, system_offset, system_cmp = [], [], []
    for v in system:
        eq_line = []
        for p in system_vars:
            if p in v['params']:
                eq_line.append(v['params'][p])
            else:
                eq_line.append(0)
        system_mat.append(eq_line)
        system_offset.append(v['offset'])
        system_cmp.append(v['cmp'])

    return Guard(np.transpose(np.array(system_vars)), np.array(system_mat), np.array(system_offset), np.array(system_cmp))
    
comparator_regex = re.compile("(<|>|=)")
guard_offset_regex = re.compile("([+-] *[0-9]+) *(<|>|=)")

def parse_guard(eq):
    '''
        Parse a single guard equation, of the form ax + by + ... + c cmp 0
    '''
    # First, find the comparator 
    m = re.search(comparator_regex, eq)
    if not m:
        raise Exception('Invalid guard syntax', eq)
    cmp = m.group(1)

    # Then determine linear parameters
    params = {}
    for mm in re.finditer(eq_params_regex, eq):
        params[mm.group(2)] = float(mm.group(1).replace(' ', ''))

    # Finally, find offset
    offset = 0
    m = re.search(guard_offset_regex, eq)
    if m:
        offset = float(m.group(1).replace(' ',''))  
    
    return {'cmp': cmp, 'params': params, 'offset': offset }


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

################################################################################

parseFile('example.json')