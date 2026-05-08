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
        if W is not None:
            x_arr = np.asarray(x)
            y_arr = np.asarray(y)
            keep = np.array(
                [W.inWindow(xi, yi) for xi, yi in zip(x_arr, y_arr)], dtype=bool
            )
            if not keep.all():
                x = list(x_arr[keep])
                y = list(y_arr[keep])
                if d is not None:
                    d = list(np.asarray(d)[keep])
                if isinstance(marks, pd.DataFrame):
                    marks = marks.iloc[keep].reset_index(drop=True)

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
