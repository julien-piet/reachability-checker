# Reachability calculation for linear first order differential equations

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from scipy.linalg import expm, sinm, cosm, det

def norme(vector):
    return np.sqrt(sum(x**2 for x in vector))
    

def lin_solver(eq1, v1, dir, values):
    '''
        returns the coordinates of the intersection point, if it exists
    '''
    N = len(eq1)
    M = np.zeros([[1 if j == i else 0 for j in range(N)] for i in range(N)])
    M[dir] = eq1
    b = np.array(values)
    b[dir] = v1
    
    try:
        return np.linalg.solve(M,b)
    except LinAlgError:
        return None
    

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
        
    def intersect_box(self, other):
        '''
            returns the intersection of the polyhedra with a box, approximated to the nearest box.
            method : calculate intersection of each hyperplan, approximate it as a box, intersect boxes.
        '''
        rtn = None
        min_p = other.min_p
        max_p = other.max_p
        for i in range(self.m):
            for p in other.points:
                
                
        


class box:
    '''
        representation of a box
    '''
    
    def __init__(self, a, b):
        self.min_p = np.array([min(a[i], b[i]) for i in range(len(a))])
        self.max_p = np.array([max(a[i], b[i]) for i in range(len(a))])
        self.points = []
        for i in range(2**len(self.min_p)):
            changes = [i & 2**j for j in range(len(self.min_p))]
            self.points.append(np.array([self.max_p[j] + changes[j]*(self.min_p[j] - self.max_p[j]) for j in range(len(self.min_p))]))
            
    def __str__(self):
        return "[" + str(self.min_p) + "," + str(self.max_p) + "]"
            
    def union(this, other):
        '''
            returns the union of this interval and the other, and any space in between
        '''
        return box(np.array([min(this.min_p[i], other.min_p[i]) for i in range(len(this.min_p))]), np.array([max(this.max_p[i], other.max_p[i]) for i in range(len(this.max_p))]))
        
    def intersection(this, other):
        '''
            returns the intersection of this interval and the other
        '''
        mins = np.array([max(this.min_p[i], other.min_p[i]) for i in range(len(this.min_p))])
        maxs = np.array([min(this.max_p[i], other.max_p[i]) for i in range(len(this.max_p))])
        return None if any(mins[i] > maxs[i] for i in range(len(mins))) else box(mins, maxs)
        
    def multiply(self, f):
        '''
            Returns a new instance of interval, given by the multiplication of its original value by f)
        '''
        return box(np.dot(f, self.min_p), np.dot(f, self.max_p))
        
    def fuzz(self, delta):
        '''
            Modifies the interval by adding a fuzzing factor
        '''
        self.min_p -= delta[0]
        self.max_p += delta[1]

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
        #beta = self.U * ((expm(self.A * self.d) - 1) / det(self.A))
        
        current_step = 0
        for i in range(n):
            self.iteration(current_step)
            current_step += 1
    
    
    def iteration(self, current_step):
        ''' 
            Calculate Reach for the next step
        '''
        d = self.d
        print(expm(self.A))
        # First, calculate the next reach
        self.reach[current_step+1] = self.reach[current_step].multiply(expm(self.d * self.A))
        #self.reach[current_step+1].fuzz(beta)
        
        # Now, update range (optimize here - Uses zonotope)
        t = current_step * d
        large_area = box.union(self.reach[current_step], self.reach[current_step+1])
        self.range.append(box(np.append([t], large_area.min_p), np.append([t+d], large_area.max_p)))
        
        
A = np.array([[0,-1],[1,0]])
U = 0
x0 = box(np.array([0,0]), np.array([1,1]))
d = 1
n = 4

crc = solver(A,U,x0,n,d)
crc.run()

print(A)
print(x0)
for b in crc.range:
    print(b)