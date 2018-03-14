import scipy

class Equation:
    def __init__(self, V, A, B):
        self.V = V
        self.A = A
        self.B = B

    def resolve(V0):
       scipy.integrate.quad(self.compute, a, b)
        
    def compute(self, t):
        pass
        
    def __str__(self):
        return "V : " + str(self.V)  + "A : " + str(self.A)  + "B : " + str(self.B) 

