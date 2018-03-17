import numpy as np
from scipy.spatial import ConvexHull
from numpy.linalg import norm
import matplotlib.pyplot as plt

def basisVect(n, i, val=1.0):
    v = np.zeros(n)
    v[i] = val
    return v

class Zonotope:

    @staticmethod
    def ball(n, c, r):
        assert c.shape[0] == n, "Invalid dimension beetween provided n and C (in ball)"
        gens = [r * basisVect(n, i) for i in range(0, n)]
        return Zonotope(n, c, gens)

    @staticmethod
    def fromIntervals(I):
        n = len(I)
        c = np.zeros(n)
        gens = []
        for i, interval in enumerate(I):
            c[i] = (interval[0] + interval[1]) * 0.5
            gens.append(basisVect(n, i, (interval[1] - interval[0]) * 0.5))
        return Zonotope(n, c, gens)

    def __init__(self, n, c=None, gens=[]):
        self.n = n 
        self.c =  np.zeros(n) if c is None else np.array(c)
        self.gens = [np.array(g) for g in gens]
        self.points = None

    def computeP(self, erA):
        assert erA.shape[0] == self.n, "Invalid dimension beetween I and erA (in P)"
        c1 = (self.c + erA.dot(self.c)) / 2
        c2 = (self.c - erA.dot(self.c)) / 2
        gens1 = [(g + erA.dot(g)) / 2 for g in self.gens]
        gens2 = [(g - erA.dot(g)) / 2 for g in self.gens]
        return Zonotope(self.n, c1, gens1 + [c2] + gens2)
        
    def __add__(self, Z):
        assert Z.n == self.n, "Invalid dimension beetween Z1 and Z2 (in add)"
        c = self.c + Z.c
        gens = self.gens + Z.gens
        return Zonotope(self.n, c, gens)

    def mult(self, L):
        assert L.shape[0] == self.n, "Invalid dimension beetween L and Z (in mul)"
        c = L.dot(self.c)
        gens = [L.dot(g) for g in self.gens]
        return Zonotope(self.n, c, gens)

    def reduce(self, m):
        while len(self.gens) > m*self.n:
            h, h2 = self.heuristicNormDiff(), []
            assert len(h) == 2*self.n, "Heuristic raised invalid reducable generator list, must be of size 2*n"

            for j in range(self.n):
                s = 0
                for i in range(2*self.n):
                    s += np.abs(self.gens[h[i]][j])
                h2.append(np.array([(s if i == j else 0) for i in range(self.n)]))

            for j in range(self.n):
                self.gens[h[j]] = h2[j]

            for j in range(self.n):
                del self.gens[h[self.n+j]]

    def heuristicNormDiff(self):
        norms = [norm(g, 1) - norm(g, np.inf) for g in self.gens]
        sortedI = np.argsort(norms)
        return sortedI[:2*self.n]

    def __str__(self):
        return "c: " + str(self.c) + "\n gens: " + str(self.gens)

    def getPoints(self):
        if not self.points is None:
            return self.points

        self.points = []
        for m in range(2**len(self.gens)):
            point = self.c
            for i, gi in enumerate(self.gens):
               point = point + gi if ((m >> i) & 1) else point - gi
            self.points.append(point)
        return self.points

    def plot(self, c):
        points = self.getPoints() 
        hull = ConvexHull(points)
        x, y = [points[i][0] for i in hull.vertices], [points[i][1] for i in hull.vertices]
        x.append(points[hull.vertices[0]][0])
        y.append(points[hull.vertices[0]][1])
        plt.plot(x, y, c=c)
        plt.draw()
        plt.pause(0.001)

    def supNorm(self):
        points = self.getPoints()
        s = 0
        for p in points:
            n = norm(p, np.inf)
            if n > s: s = n
        return s

    def is_intersecting(self, guard):
        for c in guard.C:
            assert c=='=', "Can't use zonotope on inequality guards. Only equality guards are supported"
            
        v = guard.B - guard.A.dot(self.c)
        ns, ps = np.zeros(guard.A.shape[0]), np.zeros(guard.A.shape[0])
        for g in self.gens:
            x = guard.A.dot(g)
            ps = ps + x
            ns = ns - x

        for i in range(guard.A.shape[0]):
            if x[i] > ps[i]:
                return False
            if x[i] < ns[i]:
                return False

        return True


