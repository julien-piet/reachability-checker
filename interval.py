class Interval:
    def __init__(self, a=0, b=0):
        self.a = a
        self.b = b

    def fromMax(self, m):
        if m < 0:
            m = -m
        self.a = -m
        self.b = m

    def fromInterval(self, s, a, b):
        if s == '-':
            a = -a
            b = -b
        if a > b:
            a, b = b, a
        self.a = a
        self.b = b

    def toList(self):
        return [self.a, self.b]

    def __str__(self):
        return '[' + self.a + ';' + self.b + ']'

    def __repr__(self):
        return '[' + self.a + ';' + self.b + ']'
