import scipy

class Equation:
    '''
        Default interface for an equation.
        Should always implement the fromSystem(A, B) static method returning an Equation object representing
            V' = AV + B
    '''

    @staticmethod
    def fromSystem(A, B):
        return Equation(A, B)

    def __init__(self, A, B):
        self.A = A
        self.B = B

    def __str__(self):
        return "A : " + str(self.A)  + "B : " + str(self.B) 

