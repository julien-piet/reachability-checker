# Reachability calculation for linear first order differential equations
# Diff equation doesn't work, have to review equations
### Two variables have to be declared = T (time) and X

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from scipy.spatial.qhull import QhullError
from numpy.linalg.linalg import LinAlgError

from ..generics.update     import Update
from ..generics.transition import Transition


class Equation:
    '''
        Represents an equation of the form
        X_dot = AX + [delta-, delta+]
    '''
    
    @staticmethod
    def fromSystem(A, B):
        return Equation(A[1][1], B[1][0], B[1][1])

    def __init__(self, a, d_m, d_M):
        self.a = a
        self.delta = np.array([d_m, d_M])

    def __str__(self):
        return "A : " + str(self.a)  + "Delta : " + str(self.delta) 


    
class Guard:
    '''
        Represents a simple guard
        at + bx cmp c
        OR
        AX cmp c
    '''
    
    @staticmethod
    def fromSystem(A, B, C):
        return Guard(A[0][0], A[0][1], -B[0], C[0])

    def __str__(self):
        return str(self.a) + "t + " + str(self.b) + "x " + str(self.cmp) + " " + str(self.c)
    
    def __init__(self, a, b, c, cmp):
        self.a = a
        self.b = b
        self.c = c
        self.cmp = cmp
        
        self.A = np.array([a,b])
    
    def intersect(self, zone):
        return zone.intersect_guard(self)
    
    def build_from_poly(poly):
        return guard(poly.A[0][0], poly.A[0][1], poly.B[0], "<")
        
    def is_satisfied(self, poly):
        '''
            Returns true is poly intersects with guard
        '''
        return poly.intersect_guard(self) is not None
    
    
class Polygon:
    '''
        Represents a Polygon
    '''
    
    def __init__(self, points):
        if len(points) == 0:
            self.points = None
            return
        if len(points) == 1:
            self.points = [points[0],points[0]]
            return
        try:
            convex_hull = ConvexHull(points)
            self.points = []
            for i in convex_hull.vertices:
                self.points.append(np.array(points[i]))
        except QhullError:
            points.sort(key=lambda pt: pt[0])
            max_val = max(np.array(points)[:,1])
            min_val = min(np.array(points)[:,1])
            max_point = [v for v in points if v[1] == max_val][0]
            min_point = [v for v in points if v[1] == min_val][-1]
            self.points = [np.array(max_point), np.array(min_point)]
        except ValueError:
            self.points = [points[0], points[0]]
        self.current = 0
            
    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= len(self.points):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            return self.points[self.current-1]
         
    @staticmethod
    def fromIntervals(I):
        '''
            Builds a Polygon from a list of intervals
        '''   
        min_p = np.array([i[0] for i in I])
        max_p = np.array([i[1] for i in I])
        points = []
        for i in range(2**len(I)):
            changes = [(i & 2**j) / 2**j for j in range(len(I))]
            points.append(np.array([max_p[j] + changes[j]*(min_p[j] - max_p[j]) for j in range(len(I))]))
        return Polygon(points)
        
    
    def intersect_guard(self, gd, eps=0.0001, DEBUG=False):
        '''
            Returns the intersection with a hyperplan
        '''
        if gd is None:
            return self
        inter_points = []
        for idx, val in enumerate(self.points):
            # Use matrix equation to solve intersection
            M = np.zeros((2,2))
            M[0] = gd.A
            b = np.zeros(2)
            b[0] = gd.c
            
            # Equation of the form ct + dx = e
            u = self.points[idx] - self.points[idx-1]
            c = -u[1]
            d = u[0]
            e = np.array([c,d]).dot(self.points[idx])
            
            M[1] = np.array([c,d])
            b[1] = e
            
            potential_solution = None
            try:
                potential_solution = np.linalg.solve(M,b)
            except LinAlgError:
                continue
            
            # If we get here, check if solution is valid
            max_values = [max(self.points[idx-1][i], self.points[idx][i]) for i in range(2)]
            min_values = [min(self.points[idx-1][i], self.points[idx][i]) for i in range(2)]
            if all(potential_solution[i] >= min_values[i]-eps and potential_solution[i] <= max_values[i]+eps for i in range(2)):
                inter_points.append(potential_solution)
        
        if gd.cmp == "<":
            inter_points += [value for value in self.points if np.dot(value, gd.A) <= gd.c]
        elif gd.cmp == ">":
            inter_points += [value for value in self.points if np.dot(value, gd.A) >= gd.c]
        if len(inter_points) < 2:
            return None
        
        return Polygon(inter_points)
        
    
    def __str__(self):
        return str(self.points)
            
    def union(this, other):
        '''
            returns the union of this Polygon and the other, and any space in between
        '''
        return Polygon(this.points + other.points)
        
    def adjacent(self, other, eps=0.0001):
        '''
            returns True if both polygons share at least one point, with error eps
        '''
        for i in self.points:
            for j in other.points:
                if (np.linalg.norm(i-j) < eps):
                    return True
        return False
        
    def multiply(self, f):
        '''
            Returns a new instance of interval, given by the multiplication of its original value by f)
        '''
        return Polygon([np.dot(value, f) for value in self.points])
        
    def fuzz(self, delta):
        '''
            Modifies the interval by adding a fuzzing factor
        '''
        return Polygon([value + delta[0] for value in self.points] + [value + delta[1] for value in self.points])

    def max_point(self):
        '''
            Return a point with largest second coordinate
        '''
        self.points.sort(key=lambda pt:pt[0])
        max_coord = max(v[1] for v in self.points)
        return [v for v in self.points if v[1] == max_coord][-1]
        
    def min_point(self):
        '''
            Return a point with smallest second coordinate
        '''
        self.points.sort(key=lambda pt:pt[0])
        min_coord = min(v[1] for v in self.points)
        return [v for v in self.points if v[1] == min_coord][0]
        
    def update(self, update):
        pts = []
        for v in self.points:
            pts.append(update.A @ v + update.B)
        return Polygon(pts)
        
    def plot(self, c):
        points = self.points
        t, x = [pt[0] for pt in points], [pt[1] for pt in points]
        t.append(t[0])
        x.append(x[0])
        plt.plot(t, x, c=c)
        plt.draw()
        plt.pause(0.001)
        

class solver:
    '''
        Solver main class
        Equation :
            X' = AX + U(t)
        A : scalar
        U(t) : random noize (Max abs value U)
        
        initial value in interval x0
        x0 positive
        
        For steps of duration d until guard is not satisfied OR time > T
    '''
    
    @staticmethod
    def get_set_type():
        return Polygon

    @staticmethod
    def get_equation_type():
        return Equation

    @staticmethod
    def get_guard_type():
        return Guard

    @staticmethod
    def get_update_type():
        return Update
        
    def __init__(self, d, debug=True):
        self.debug = debug
        self.d = d
        
    def set_params(self, node, T, init_val, color):
        # User defined variables
        self.a = node.equation.a
        self.eq = node.equation
        self.guard = node.guard
        self.links = node.links
        self.T = T
        
        # Derived variables
        self.n = int(T / self.d)
        self.reach = [0 for i in range(self.n+1)]
        self.reach[0] = init_val
        self.range = []
        
        eq = node.equation
        self.beta = np.array([[self.d, eq.delta[0] * ((np.exp(eq.a * self.d) - 1) / eq.a)],[self.d, eq.delta[1] * ((np.exp(eq.a * self.d) - 1) / eq.a)]])
        
        self.reach[0].plot(c="red")
        
    
    def solve(self, eq, init_val, node, T, color):
        
        # User defined variables
        self.set_params(node, T, init_val, color)
        current_step = 0
        trs = []
        for i in range(self.n):
            self.iteration(current_step, self.beta)
            current_step += 1
            
            # Intersecting with node guard
            self.range[-1] = self.range[-1].intersect_guard(self.guard)
            if self.range[-1] is None:
                self.range.pop()
                break
            
            # Finding new nodes
            self.range[-1].plot(color)
            for l in node.links:
                new_vals = self.range[-1].intersect_guard(l.guard)
                if new_vals:
                    trs.append(Transition(l, new_vals.update(l.update)))
        trs = self.simplify(trs)
        return trs
    
    
    def simplify(self, trs):
        '''
            Simplifies transitions to only keep connex parts (over-approximates)
        '''
        # First, build connext components (We could use union-find for better efficiency)
        elts = []
        for elt in trs:
            saved = False
            for set in elts:
                for val in set:
                    if val.values.adjacent(elt.values):
                        set.append(elt)
                        saved = True
                        break
                if saved:
                    break
            if not saved:
                elts.append([elt])
        
        # Now, union
        trs = []
        for connex_set in elts:
            new_set = connex_set[0]
            for set in connex_set:
                new_set.values = Polygon.union(new_set.values,set.values)
            trs.append(new_set)
        return trs
   
    def iteration(self, current_step, beta):
        ''' 
            Calculate Reach for the next step
        '''
        d = self.d
        # First, calculate the next reach
        self.reach[current_step+1] = self.reach[current_step].multiply(np.array([[1, 0],[0, np.exp(self.d * self.a)]])).fuzz(beta)
        
        # Now, update range (optimize here - Uses polygons)
        # Polygon calculation
        # secant --- y : gamma * t + b
        # tangent -- y : gamma * t + phi
        min_point = self.reach[current_step].min_point()
        start = min_point[1]
        A = self.a
        t = min_point[0]
        gamma = start * (np.exp(A*d) - 1)  / d
        b = start - (gamma * t)
        phi = gamma * (1 - A*t - np.log((gamma / (start * A)))) / A if start != 0 else b

        #Adding to range
        pts = []
        max_point = self.reach[current_step].max_point()
        pts.append(max_point)
        pts.append([t, phi + (gamma * t) + beta[1][1]])
        pts.append([t+self.d, phi + (gamma * (t+self.d)) + beta[1][1]])
        pts.append(self.reach[current_step+1].max_point())
        self.range.append(Polygon(pts))