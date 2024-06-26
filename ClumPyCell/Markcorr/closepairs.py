import logging
import math

import numpy as np


def closepairs(
    X,
    rmax,
    d=None,
    twice=True,
    what=None,
    distinct=True,
    neat=True,
    periodic=True,
    pp=None,
):
    if rmax < 0:
        raise ValueError("rmax must be non-negative")
    if X.getWindow().getType() != "rectangular":
        logging.warning("Periodic edge correction applied in non-rectangular window")

    npts = len(X.getX())
    if not what:
        what = "all"
    null_answer = {
        "i": [],
        "j": [],
        "xi": [],
        "yi": [],
        "xj": [],
        "yj": [],
        "dx": [],
        "dy": [],
        "d": [],
        "Awt": [],
    }

    if npts == 0:
        return null_answer

    xsort = np.sort(X.getX())
    ysort = np.sort(X.getY())
    npairs = npts * npts

    if npairs <= 1024:
        nsize = 1024
    else:
        catchfraction = math.pi * rmax**2 / X.getWindow().getArea()
        nsize = min(max(1024, int(2 * catchfraction * npairs)), npairs)

    if periodic:
        x = X.getX()
        y = X.getY()
        d = X.getD() if d is not None else np.zeros(len(x))
        z = closePpairs(x, y, d, rmax, nsize, pp)
        i, j, d, areaWt = z[0], z[1], z[2], z[3]

        if what == "all":
            xi, yi = x[np.array(i) - 1], y[np.array(i) - 1]
            xj, yj = x[np.array(j) - 1], y[np.array(j) - 1]
            dx, dy = xj - xi, yj - yi
    else:
        if not distinct:
            null_answer = {
                "i": list(range(npts)),
                "j": list(range(npts)),
                "xi": X.getX(),
                "yi": X.getY(),
                "xj": X.getX(),
                "yj": X.getY(),
                "dx": [0] * npts,
                "dy": [0] * npts,
                "d": [0] * npts,
            }

        z = Fclosepairs(npts, xsort, ysort, rmax)
        npairs = z[-1]

        if npairs <= 0:
            return null_answer

        i, j = z[0][:npairs], z[1][:npairs]
        if what == "all":
            xi, yi = z[2][:npairs], z[4][:npairs]
            xj, yj = z[3][:npairs], z[5][:npairs]
            dx, dy = z[6][:npairs], z[7][:npairs]
            d = z[8][:npairs]

    if twice:
        i, j = np.concatenate([i, j]), np.concatenate([j, i])
        if what == "all":
            xi, yi = np.concatenate([xi, xj]), np.concatenate([yi, yj])
            xj, yj = np.concatenate([xj, xi]), np.concatenate([yj, yi])
            dx, dy = np.concatenate([dx, -dx]), np.concatenate([dy, -dy])
            d = np.concatenate([d, d])
    else:
        if neat:
            swap = np.array(i) > np.array(j)
            i[swap], j[swap] = j[swap], i[swap]
            if what == "all":
                xi[swap], xj[swap] = xj[swap], xi[swap]
                yi[swap], yj[swap] = yj[swap], yi[swap]
                dx[swap], dy[swap] = -dx[swap], -dy[swap]

    if what == "all":
        answer = {
            "i": i,
            "j": j,
            "xi": xi,
            "yi": yi,
            "xj": xj,
            "yj": yj,
            "dx": dx,
            "dy": dy,
            "d": d,
            "Awt": areaWt,
        }
    elif what == "indices":
        answer = {"i": i, "j": j}
    elif what == "ijd":
        answer = {"i": i, "j": j, "d": d}

    return answer


def paircount(nxy, x, y, rmaxi):
    count = 0
    r2max = rmaxi * rmaxi
    if nxy == 0:
        return 0
    maxchunk = 0
    while maxchunk < nxy:
        for i in range(maxchunk, min(maxchunk + 65536, nxy)):
            xi, yi = x[i], y[i]
            for j in range(i + 1, nxy):
                dx, dy = x[j] - xi, y[j] - yi
                if dx * dx + dy * dy <= r2max:
                    count += 1
        maxchunk += 65536
    return count


def closePpairs(xx, yy, dia, rr, nguess, pp=None):
    n = len(xx)
    r2max = rr * rr
    pp_x = pp.getX() if pp else []
    pp_y = pp.getY() if pp else []
    pp_d = pp.getD() if pp else []
    jout, iout, dout, areaWt = [], [], [], []

    if n > 0 and nguess > 0:
        for i in range(n):
            xi, yi, di = xx[i], yy[i], dia[i]
            for j in range(i + 1, n):
                dx, dy = abs(xx[j] - xi), abs(yy[j] - yi)
                if dx < rr and dy < rr:
                    d2 = max(0, math.sqrt(dx * dx + dy * dy) - di / 2 - dia[j] / 2)
                    distance = math.sqrt(d2)
                    area = overlapA(di, distance, distance)
                    for k in range(len(pp_x)):
                        p1, p2, p3 = (
                            np.array([xi, yi]),
                            np.array([xx[j], yy[j]]),
                            np.array([pp_x[k], pp_y[k]]),
                        )
                        d = np.linalg.norm(np.cross(p2 - p1, p1 - p3)) / np.linalg.norm(
                            p2 - p1
                        )
                        if d < pp_d[k] / 2:
                            dist_pass = 2 * math.sqrt(pp_d[k] ** 2 / 4 - d**2)
                            if dist_pass < d2:
                                d2 -= dist_pass
                    if d2 * d2 <= r2max:
                        jout.append(j + 1)
                        iout.append(i + 1)
                        dout.append(d2)
                        areaWt.append(area)
    return iout, jout, dout, areaWt


def overlapA(r1, r2, d):
    if d >= r1 + r2:
        return 0
    elif d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    alpha1 = math.acos((d**2 + r1**2 - r2**2) / (2 * d * r1))
    alpha2 = math.acos((d**2 + r2**2 - r1**2) / (2 * d * r2))
    A1 = r1**2 * alpha1
    A2 = r2**2 * alpha2
    A3 = -0.5 * math.sqrt(
        (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)
    )
    return A1 + A2 + A3


def Fclosepairs(nxy, x, y, r):
    r2max = r * r
    jout, iout, xjout, xiout, yjout, yiout, dxout, dyout, dout = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    maxchunk = 0
    while maxchunk < nxy:
        for i in range(maxchunk, min(maxchunk + 65536, nxy)):
            xi, yi = x[i], y[i]
            for j in range(i + 1, nxy):
                dx, dy = x[j] - xi, y[j] - yi
                if dx * dx + dy * dy <= r2max:
                    jout.append(j + 1)
                    iout.append(i + 1)
                    xjout.append(x[j])
                    xiout.append(xi)
                    yjout.append(y[j])
                    yiout.append(yi)
                    dxout.append(dx)
                    dyout.append(dy)
                    dout.append(math.sqrt(dx * dx + dy * dy))
        maxchunk += 65536
    return iout, jout, xiout, xjout, yiout, yjout, dxout, dyout, dout, len(jout)
