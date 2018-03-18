class UnionFind:

    def union(p, q):
        i = UnionFind.find(p)
        j = UnionFind.find(q)
        if i != j:
            q.parent = i
        
    def find(p):
        if p.parent is not p:
            p.parent = UnionFind.find(p.parent)
        return p.parent
        

class Elt:
    
    def __init__(self, value):
        self.value = value
        self.parent = self