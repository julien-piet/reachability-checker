# Reachability calculation for linear first order differential equations

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from scipy.spatial.qhull import QhullError
class interval:
    '''
        representation of an interval
        [a, b]
    '''
    
    def __init__(self, a, b):
        if a <= b:
            self.a = a
            self.b = b
            self.void = False
        else:
            self.void = True
            
    def __str__(self):
        return "[" + str(self.a) + "," + str(self.b) + "]"
            
    def union(this, other):
        '''
            returns the union of this interval and the other, and any space in between
        '''
        if other.void:
            return interval(this.a, this.b)
        if this.void:
            return interval(other.a, other.b)
        
        return interval(min(this.a, other.a), max(this.b, other.b))
    
    def union_set(intervals):
        '''
            Returns the union of all given intervals. Non-destructive
        '''
        union = interval(intervals[0].a, intervals[0].b)
        for elt in intervals[1:]:
            union = interval.union(elt, union)
        return union
        
    def multiply(self, f):
        '''
            Returns a new instance of interval, given by the multiplication of its original value by f)
        '''
        if f > 0:
            return interval(self.a * f, self.b * f)
        elif f < 0:
            return interval(self.b * f, self.a * f)
        return interval(0,0)
        
    def fuzz(self, delta):
        '''
            Modifies the interval by adding a fuzzing factor
        '''
        self.a = self.a - delta
        self.b = self.b + delta

class guard:
    '''
        Represents a simple guard
        at + bx cmp c
        OR
        AX cmp c
    '''
    
    def __init__(self, a, b, c, cmp):
        self.a = a
        self.b = b
        self.c = c
        self.cmp = cmp
        
        self.A = np.array([a,b])
    
    def build_from_poly(poly):
        return guard(poly.A[0][0], poly.A[0][1], poly.B[0], "<")
    
    
class polygon:
    '''
        Represents a polygon
    '''
    
    def __init__(self, points):
        try:
            convex_hull = ConvexHull(points)
            self.points = []
            for i in convex_hull.vertices:
                self.points.append(np.array(points[i]))
        except QhullError:
            max_val = max(np.array(points)[:,1])
            min_val = min(np.array(points)[:,1])
            max_point = [v for v in points if v[1] == max_val][0]
            min_point = [v for v in points if v[1] == min_val][0]
            self.points = [np.array(max_point), np.array(min_point)]
            
    
    def intersect_guard(self, guard):
        '''
            Returns the intersection with a hyperplan
        '''
        inter_points = []
        for idx, val in enumerate(self.points):
            # Use matrix equation to solve intersection
            M = np.zeros((2,2))
            M[0] = guard.A
            b = np.zeros(2)
            b[0] = guard.B
            
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
            if potential_solution[dir] >= min_values and  potential_solution[dir] <= max_values:
                inter_points.append(potential_solution)
        if len(inter_points) == 0:
            return None
        
        if guard.cmp == "<":
            inter_points += [value for value in self.points if np.dot(value, guard.A) <= guard.c]
        elif guard.cmp == ">":
            inter_points += [value for value in self.points if np.dot(value, guard.A) >= guard.c]
        return polygon(inter_points)
        
    
    def __str__(self):
        return "STR method not written"
            
    def union(this, other):
        '''
            returns the union of this polygon and the other, and any space in between
        '''
        return polygon(this.points + other.points)
        
    def multiply(self, f):
        '''
            Returns a new instance of interval, given by the multiplication of its original value by f)
        '''
        return polygon([np.dot(value, f) for value in self.points])
        
    def fuzz(self, delta):
        '''
            Modifies the interval by adding a fuzzing factor
        '''
        return polygon([value + delta[0] for value in self.points] + [value + delta[1] for value in self.points])

    def max_point(self):
        '''
            Return a point with largest second coordinate
        '''
        max_coord = max(v[1] for v in self.points)
        return [v for v in self.points if v[1] == max_coord][0]
        
    def min_point(self):
        '''
            Return a point with smallest second coordinate
        '''
        min_coord = min(v[1] for v in self.points)
        return [v for v in self.points if v[1] == min_coord][0]

class solver:
    '''
        Solver main class
        Equation :
            X' = AX + U(t)
        A : scalar
        U(t) : random noize (Max abs value U)
        
        initial value in interval x0
        x0 positive (TBF)
        
        For n steps of duration d
    '''
    
    def __init__(self, A, U, x0, n, d):
        # User defined variables
        
        self.A = A
        self.U = U
        self.n = n
        self.d = d
        
        # Derived variables
        
        self.T = d*n
        self.reach = [0 for i in range(n+1)]
        self.reach[0] = x0
        
        #range will contain points of a convex enveloppe of solutions
        self.range = []
        
    
    def run(self):
        '''
            Run the algorithm and store reachabilities in array self.reach
        '''        
        #Calculate beta, the bloating factor
        beta = np.array([[self.d, self.U * ((np.exp(self.A * self.d) - 1) / self.A)],[self.d, -self.U * ((np.exp(self.A * self.d) - 1) / self.A)]])

        current_step = 0
        for i in range(n):
            self.iteration(current_step, beta)
            current_step += 1
    
    
    def iteration(self, current_step, beta):
        ''' 
            Calculate Reach for the next step
        '''
        d = self.d
        # First, calculate the next reach
        self.reach[current_step+1] = self.reach[current_step].multiply(np.array([[1, 0],[0, np.exp(self.d * self.A)]])).fuzz(beta)
        
        # Now, update range (optimize here - Uses zonotope)
        # Zonotope calculation
        # secant --- y : gamma * t + b
        # tangent -- y : gamma * t + phi
        min_point = self.reach[current_step].min_point()
        start = min_point[1]
        A = self.A
        t = current_step * d
        gamma = start * (np.exp(A*d) - 1)  / d
        b = start - (gamma * t)
        phi = gamma * (1 - A*t - np.log((gamma / (start * A)))) / A
    
        #Adding to range
        pts = []
        max_point = self.reach[current_step].max_point()
        pts.append(max_point)
        pts.append([min_point[0], phi + (gamma * t) + beta[1][1]])
        pts.append([self.reach[current_step+1].min_point()[0], phi + (gamma * (t+delta)) + beta[1][1]])
        pts.append(self.reach[current_step+1].max_point())
        self.range.append(polygon(pts))
        
        
a = 1
u = 0.1
x0 = polygon([[0,2],[0,3]])
delta = 0.1
n = 100
slv = solver(a,u,x0,n,delta)
slv.run()
print(slv.range)
# Visualisation
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_xlim([0, 4])
ax.set_ylim([0, 3*np.exp(4)])
for i in slv.range:
    polygon= plt.Polygon(i.points, edgecolor='r')
    ax.add_patch(polygon)
    
# Exponential for comparison
time = np.linspace(0,4,50)
plt.plot(time, slv.reach[0].min_point()[1]*np.exp(a*time), color='y')
plt.plot(time, slv.reach[0].max_point()[1]*np.exp(a*time), color='y')
plt.show() 
