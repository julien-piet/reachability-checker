import scipy

class Equation:
    @staticmethod
    def fromSystem(A, B):
        return Equation(A, B)

    def __init__(self, A, B):
        self.A = A
        self.mu = 0
        for l in B:
            for c in l:
                if abs(c) > self.mu:
                    self.mu = abs(c)

    def __str__(self):
        return "A : " + str(self.A)  + "mu : " + str(self.mu) 

