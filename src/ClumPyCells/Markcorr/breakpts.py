import logging, sys, math
import statistics
import numpy as np


class breakpts:
    def __init__(self, val, maxi, even=False, npos=None, step=None):
        self.val = val
        self.maxi = maxi
        self.even = even
        self.npos = npos
        self.step = step
        self.ncells = len(val) - 1

    def getR(self):
        return self.val

    def getMax(self):
        return self.maxi


def handle_r_b(r=None, breaks=None, window=None, pixeps=None, rmaxdefault=None):
    if r and breaks:
        logging.error("Do not specify both " + r + " and " + breaks)
        sys.exit(2)
    if breaks:
        return breaks
    elif r:
        return breakpts_from_r(r)
    else:
        pixeps = rmaxdefault / 128
        if not pixeps:
            if window.getType() == "mask":
                pixeps = window.minEdge()
            else:
                pixeps = rmaxdefault / 128
        rstep = pixeps / 4
        breaks = make_even_breaks(rmaxdefault, rstep)
        return breaks


def breakpts_from_r(r):
    if not r:
        logging.error("Input r is not specified")
        sys.exit(2)
    if len(r) < 2:
        logging.error("r has length", len(r), "- must be at least 2")
        sys.exit(2)
    if r[0] != 0:
        logging.error("First r value must be 0")
        sys.exit(2)
    res = all(i < j for i, j in zip(r, r[1:]))
    if not res:
        logging.error("successive values of r must be increasing")
        sys.exit(2)
    dr = r[1] - r[0]
    r.insert(0, -dr)
    return r


def make_even_breaks(bmax, bstep, npos=None):
    if bmax <= 0:
        logging.error("bmax must be positive")
        sys.exit(2)
    if npos:
        bstep = bmax / npos
        val = list(range(0, bmax, bmax / npos))
        val.insert(0, -bstep)
        right = bmax
    else:
        npos = math.ceil(bmax / bstep)
        right = bstep * npos
        val = list(np.arange(0, float(right), bstep))
        val.insert(0, -bstep)
    return breakpts(val, right, True, npos, bstep)


def as_breakpts(b):
    if len(b) > 2:
        if b[1] != 0:
            logging.error("breakpoints do not satisfy breaks[1] = 0")
            sys.exit(2)
        steps = [b[i + 1] - b[i] for i in range(len(b) - 1)]
        if max(steps) - min(steps) < 1e-07 * statistics.mean(steps):
            return breakpts(b, max(b), True, len(b) - 2, steps[0])
        else:
            return breakpts(b, max(b), False)
    return
