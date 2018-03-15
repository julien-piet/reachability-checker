import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from numpy.linalg import norm

from ..generics.equation import Equation
from ..generics.guard    import Guard
from ..generics.update   import Update
from .zonotope           import Zonotope

class Solver:

    @staticmethod
    def get_set_type():
        return Zonotope

    @staticmethod
    def get_equation_type():
        return Equation

    @staticmethod
    def get_guard_type():
        return Guard

    @staticmethod
    def get_update_type():
        return Update


    def __init__(self, r, m, debug=True):
        self.debug = debug
        self.r = r
        self.m = m

    def solve(self, eq, init, guards):
        pass

    def set_params(self, A, I, mu):
        self.A = np.array(A)
        self.n = self.A.shape[0]
        self.I = I
        self.mu = mu
        if self.debug:
            print("---------------------")
            print("A =\n", self.A)
            print("I =\n", self.I)
            print("r =", self.r)
            print("mu =", self.mu)
            print("m =", self.m)
        
        self.erA = expm(self.A * self.r)
        self.nA = norm(self.A, np.inf)
        if self.debug:
            print("---------------------")
            print("exp(rA) =\n", self.erA)
            print("||A|| =", self.nA)
        
        self.betaR = (np.exp(self.r * self.nA)-1)/self.nA * self.mu
        self.alphaR = (np.exp(self.r * self.nA) - 1 - self.r*self.nA) * self.I.supNorm()
        if self.debug:
            print("---------------------")
            print("betaR =", self.betaR)
            print("alphaR =", self.alphaR)

        self.abBall = Zonotope.ball(self.n, np.zeros(self.n), self.alphaR + self.betaR)
        self.bBall  = Zonotope.ball(self.n, np.zeros(self.n), self.betaR)
        if self.debug:
            print("---------------------")
            print("ball(alphaR + betaR) =\n", self.abBall)
            print("ball(betaR) =\n", self.bBall)

        self.P = I.computeP(self.erA)
        self.Q = self.P + self.abBall
        self.stepno = 0
        self.t = 0
        if self.debug:
            print("---------------------")
            print("Step", self.stepno, "(t =", self.t, ")")
            print("P =\n", self.P)
            print("Q =\n", self.Q)


    def step(self):
        self.stepno += 1
        self.t += self.r
        self.P = self.Q.mult(self.erA)
        self.Q = self.P + self.bBall
        self.Q.reduce(self.m)
        if self.debug:
            print("---------------------")
            print("Step", self.stepno, "(t =", self.t, ")")
            print("P =\n", self.P)
            print("Q =\n", self.Q)



# I = Zonotope(2, [1.0, 0.0], [[0.1, 0.0], [0.0, -0.1]])
# A = [[-1, -4], [4, -1]]

# RZ = ReachabilityZono(A, I, 0.02, 0.05, 7)

# for i in range(100):
    # RZ.Q.plot(plt)
    # RZ.step()

# plt.show()
