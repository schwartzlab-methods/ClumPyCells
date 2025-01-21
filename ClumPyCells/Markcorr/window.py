class window ():
    def __init__(self, xrange, yrange, type = "rectangular", bdry = None) -> None:
        if type == "rectangular":
            self.xrange = xrange
            self.yrange = yrange
        self.type = type

    def inWindow(self, x, y):
        if x >= self.xrange[0] and x <= self.xrange[1] and y >= self.yrange[0] and y <= self.yrange[1]:
            return True
        else:
            return False
    
    def getXrange(self):
        return self.xrange

    def getYrange(self):
        return self.yrange

    def minEdge(self):
        return min(self.xrange[1]-self.xrange[0],self.yrange[1]-self.yrange[0])

    def getType(self):
        return self.type
    
    def getArea(self):
        if self.type == "rectangular":
            return max(self.xrange)*max(self.yrange)

    def getSideLength(self):
        if self.type == "rectangular":
            return [self.xrange[1]-self.xrange[0],self.yrange[1]-self.yrange[0]]

def isSubWindow(A, B):
    if A.xrange[0] >= B.xrange[0] and A.xrange[1] <= B.xrange[1] and A.yrange[0] >= B.yrange[0] and A.yrange[1] <= B.yrange[1]:
        return True
    else:
        return False
