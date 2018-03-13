# Reachability calculation for linear first order differential equations

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


def norme(vector):
    return np.sqrt(sum(x**2 for x in vector))

class polyhedra:
    '''
        representation of a n dimentional polygon with its guards
        P = { x in Rn | Ax <= B }
        A in R(m*n), B in Rm
    '''
    
    def __init__(self, m=0, n=1, A=[[]], B=[[]]):
        self.A = A
        self.B = B
        self.m = m
        self.n = n
        
    def build_from_vect(vector_list, guard_list):
        if len(vector_list) != len(guard_list) or len(vector_list) == 0:
            return None
        m = len(guard_list)
        n = len(vector_list[0])
        A = np.zeros((m,n))
        B = np.zeros(m)
        
        for i in range(m):
            A[i] = vector_list[i]
            B[i] = guard_list[i]
        
        return polyhedra(m,n,A,B)
        
    def build_from_points(point_list):
        if len(point_list) < 2:
            return None
            
        if len(point_list) == 2:
            n = len(point_list[0])
            m = 2*n
            A = np.zeros((m,n))
            B = np.zeros(m)
            #On commence par coder les bornes du segment
            point_a = point_list[0]
            point_b = point_list[1]
            p_vect = point_a - point_b
            A[0] = p_vect
            B[0] = -np.dot(p_vect,point_a)
            A[1] = -p_vect
            B[1] = np.dot(p_vect,point_b)
            current_index = 2
            
            #Ensuite, on va coder les bornes sur les normales.
            i = min([i for i,v in enumerate(p_vect) if v != 0])
            for j in range(n):
                if j == i:
                    continue
                n_vect = np.zeros(n)
                n_vect[j] = 1
                n_vect[i] = - (p_vect[i] / p_vect[j]) if p_vect[j] else 0
                A[current_index] = n_vect
                B[current_index] = np.dot(n_vect,point_a)
                A[current_index + 1] = -n_vect
                B[current_index + 1] = np.dot(n_vect,point_a)
                current_index += 2
            
            return polyhedra(m,n,A,B)
            
        convex_hull = ConvexHull(point_list)
        m = convex_hull.simplices.shape[0]
        n = convex_hull.simplices.shape[1]
        A = np.zeros((m,n))
        B = np.zeros(m)
        
        for index, equation in enumerate(convex_hull.equations):
            norm = equation[:3]
            b = -equation[3]
            A[index] = norm
            B[index] = b
        
        return polyhedra(m,n,A,B)
        
        
    def intersect_similar(self, other):
        '''
            returns a new polyhedra, formed by the intersection of both. Only works on polyhedra with same directions
        '''
        constraints = []
        values = []
        used = []
        for index,equation in enumerate(self.A):
            if all(np.linalg.det(equation,vect) != 0. for vect in other.A):
                contraints.append(equation)
                values.append(self.B[index])
            else:
                other_index = [i for i,v in other.A if np.linalg.deg(equation,v) == 0.][0]
                used.append(other_index)
                if (self.B[index] / norme(equation)) < (other.B[other_index] / norme(other.A[other_index])):
                    constraints.append(equation)
                    values.append(self.B[index])
                else:
                    constraints.append(other.A[other_index])
                    values.append(other.B[other_index])
        for index,equation in enumerate(other.A):
            if index in used:
                continue
            contraints.append(equation)
            values.append(other.B[index])
        
        A = np.array(constraints)
        B = np.array(values)
        n = A.shape[1]
        m = A.shape[0]
        return polyhedra(m,n,A,B) 
        
    


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
        beta = self.U * ((np.exp(self.A * self.d) - 1) / self.A)
        
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
        self.reach[current_step+1] = self.reach[current_step].multiply(np.exp(self.d * self.A))
        self.reach[current_step+1].fuzz(beta)
        
        # Now, update range (optimize here - Uses zonotope)
        # Zonotope calculation
        # secant --- y : gamma * t + b
        # tangent -- y : gamma * t + phi
        start = self.reach[current_step].a
        A = self.A
        t = current_step * d
        gamma = start * (np.exp(A*d) - 1)  / d
        b = start - (gamma * t)
        phi = gamma * (1 - A*t - np.log((gamma / (start * A)))) / A
        
        #Adding to range
        self.range.append([])
        self.range[-1].append([t, self.reach[current_step].b])
        self.range[-1].append([t, phi + (gamma * t) - beta])
        self.range[-1].append([t + d, phi + (gamma * (t+delta)) - beta])
        self.range[-1].append([t + d, self.reach[current_step+1].b])
        
        
a = 1
u = 0
x0 = interval(2,3)
delta = 1
n = 4
slv = solver(a,u,x0,n,delta)
slv.run()
print(slv.range)
# Visualisation
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_xlim([0, 4])
ax.set_ylim([0, 3*np.exp(4)])
for i in slv.range:
    polygon= plt.Polygon(i, edgecolor='r')
    ax.add_patch(polygon)
    
# Exponential for comparison
time = np.linspace(0,4,50)
plt.plot(time, slv.reach[0].a*np.exp(a*time), color='y')
plt.plot(time, slv.reach[0].b*np.exp(a*time), color='y')
plt.show() 
