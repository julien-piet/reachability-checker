class Guard:
    def __init__(self, vars, A, B, cmp):
        self.vars = vars
        self.cmp = cmp
        self.A = A
        self.B = B
        
    def __str__(self):
        return " Vars : " + str(self.vars) + "\n" + str(self.A) + " X + " + str(self.B) + " " + str(self.cmp) + " 0\n"