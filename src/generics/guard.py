class Guard:
    '''
        Basic interface for a Guard class.
        Should always implement a fromSystem(A, B, comp_op) static method returning a Guard object representing
            AV + B <C> 0
        Where C is a vector storing comparison operator for each row (=|<=|<|>=|>)
    '''

    @staticmethod
    def fromSystem(A, B, C):
        return Guard(A, B, C)

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C
        
    def __str__(self):
        return str(self.A) + " X + " + str(self.B) + " " + str(self.C) + " 0\n"
