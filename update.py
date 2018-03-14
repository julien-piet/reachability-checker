class Update:
    def __init__(self, vars, A, B):
        self.vars = vars
        self.A = A
        self.B = B
        
    def update(self, values):
        for v in values:
            v = A @ v + B
        return values
        
    def __str__(self):
        return " Vars : " + str(self.vars) + "\n X = " + str(self.A) + " X + " + str(self.B) + "\n"
