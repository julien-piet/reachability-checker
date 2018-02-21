# Reachability calculation for linear first order differential equations

import numpy as np
import matplotlib.pyplot as plt

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
        Amp = self.reach[0].a
        A = self.A
        t = current_step * d
        gamma = Amp * np.exp(A*t) * (np.exp(A*d) - 1)  / d
        b = (Amp * np.exp(A*t)) - (gamma * t)
        phi = gamma * (1 - np.log((gamma / (Amp * A)))) / A
        
        #Adding to range
        self.range.append([])
        self.range[-1].append([t, self.reach[current_step].b])
        self.range[-1].append([t, phi + (gamma * t) - u])
        self.range[-1].append([t + d, phi + (gamma * (t+delta)) - u])
        self.range[-1].append([t + d, self.reach[current_step+1].b])
        
        
a = 1
u = 0.1
x0 = interval(0.5,1)
delta = 1
n = 10
slv = solver(a,u,x0,n,delta)
slv.run()

# Visualisation
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_xlim([0, 10])
ax.set_ylim([0, np.exp(10)])
for i in slv.range:
    polygon= plt.Polygon(i, edgecolor='r')
    ax.add_patch(polygon)
    
# Exponential for comparison
time = np.linspace(0,10,50)
plt.plot(time, slv.reach[0].a*np.exp(a*time), color='y')
plt.plot(time, slv.reach[0].b*np.exp(a*time), color='y')
plt.show() 
