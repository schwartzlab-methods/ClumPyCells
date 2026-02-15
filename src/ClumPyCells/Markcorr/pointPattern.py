from .window import *
import pandas as pd
import logging
import numpy as np


class pointPattern:
    def __init__(self, x, y, d=None, W=None, marks=None):
        assert len(x) == len(y), "x and y coordinates are with different length"
        if isinstance(marks, pd.DataFrame):
            if len(marks) != len(x):
                logging.error("marks are with different length with coordinates")
        if W != None:
            assert isinstance(W, window)

        # filter and leave only the points within the window
        if W != None:
            for i in range(0, len(x)):
                if not W.inWindow(x[i], y[i]):
                    x.pop(i)
                    y.pop(i)
                    if isinstance(marks, pd.DataFrame):
                        marks = marks.drop(i, axis=0, inplace=True)

        self.x = np.array(x)
        self.y = np.array(y)
        self.n = len(x)
        self.marks = marks
        self.window = W
        if d is None:
            d = [0] * self.n

        self.diameter = np.array(d)

    def getMarks(self):
        return self.marks

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getD(self):
        return self.diameter

    def getWindow(self):
        return self.window

    # def getArea(self):
    #     # calculate area of the window
    #     return self.window.getArea()

    def getSideLength(self):
        return self.window.getSideLength()


def ppsubset(X, Y):
    # X is a pointPattern obj
    # Y is a list of marks
    # return a pointPattern which is part of Y belongs to the subset of X
    marx = X.getMarks()
    outputY = []
    for y in Y:
        if y in marx:
            outputY.append(y)
    return outputY
