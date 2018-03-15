class Update:
    '''
        Basic interface for an Update class.
        Should always implement a fromSystem(A, B) static method returning an Update object representing
            V = AV + B
    '''

    @staticmethod
    def fromSystem(A, B):
        return Update(A, B)

    def __init__(self, A, B):
        self.A = A
        self.B = B
        
    def update(self, values):
        for v in values:
            v = A @ v + B
        return values
        
    def __str__(self):
        return "X = " + str(self.A) + " X + " + str(self.B) + "\n"
