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
    """Enumerate close pairs of points with cell-size correction.

    The reported pairwise distance is

        d2 = max(0, ||p_i - p_j|| - r_i - r_j - sum_k chord_k)

    where ``r_i = dia_i / 2`` is the radius of point *i* (the point cells), and
    the optional ``chord_k`` term accounts for *occluder* cells (``pp``)
    intersecting the segment between ``p_i`` and ``p_j``: each occluder of
    diameter ``D_k`` whose centre projects onto the segment with perpendicular
    distance ``perp_k < D_k / 2`` removes a chord of length
    ``2 * sqrt((D_k/2)^2 - perp_k^2)`` from the corrected distance.

    The implementation is fully vectorised (NumPy) so it scales to thousands of
    points without falling back to nested Python loops.
    """
    n = len(xx)
    r2max = rr * rr
    xx = np.asarray(xx, dtype=float)
    yy = np.asarray(yy, dtype=float)
    dia = np.asarray(dia, dtype=float)

    has_pp = pp is not None and len(pp.getX()) > 0
    if has_pp:
        pp_x = np.asarray(pp.getX(), dtype=float)
        pp_y = np.asarray(pp.getY(), dtype=float)
        pp_d = np.asarray(pp.getD(), dtype=float)

    iout, jout, dout, areaWt = [], [], [], []

    if n <= 1 or nguess <= 0:
        return iout, jout, dout, areaWt

    # 1) Bounding-box pre-filter on |dx|, |dy| < rr (matches the original
    #    O(n^2) double loop semantics, including the *strict* inequality).
    i_idx, j_idx = np.triu_indices(n, k=1)
    dx = xx[j_idx] - xx[i_idx]
    dy = yy[j_idx] - yy[i_idx]
    box = (np.abs(dx) < rr) & (np.abs(dy) < rr)
    if not box.any():
        return iout, jout, dout, areaWt
    i_idx = i_idx[box]
    j_idx = j_idx[box]
    dx = dx[box]
    dy = dy[box]

    # 2) Centre-to-centre distance, then deduct the radii of the two cells.
    dist = np.sqrt(dx * dx + dy * dy)
    d2 = np.maximum(0.0, dist - dia[i_idx] / 2.0 - dia[j_idx] / 2.0)

    # 3) Subtract the chord swept inside every large occluder cell that the
    #    segment passes through. Only occluders whose foot-of-perpendicular
    #    actually lies on the segment count; otherwise the line meets the
    #    occluder outside the segment and no distance is occluded.
    if has_pp:
        p1x = xx[i_idx]
        p1y = yy[i_idx]
        seg_dx = dx
        seg_dy = dy
        seg_len2 = seg_dx * seg_dx + seg_dy * seg_dy
        # Iterate over occluders (typically much smaller than the pair count);
        # this preserves the original sequential semantics where each chord is
        # subtracted from the *running* corrected distance.
        for k in range(pp_x.size):
            r_pp = pp_d[k] / 2.0
            if r_pp <= 0.0:
                continue
            ex = pp_x[k] - p1x
            ey = pp_y[k] - p1y
            # Parameter of the foot of perpendicular along the segment, in [0,1]
            # when it lies on the segment.
            t = (ex * seg_dx + ey * seg_dy) / seg_len2
            on_seg = (t > 0.0) & (t < 1.0)
            # Perpendicular distance from occluder centre to the infinite line.
            cross = seg_dx * ey - seg_dy * ex
            perp2 = (cross * cross) / seg_len2
            intersects = on_seg & (perp2 < r_pp * r_pp)
            if not intersects.any():
                continue
            chord = np.zeros_like(d2)
            chord[intersects] = 2.0 * np.sqrt(
                np.maximum(0.0, r_pp * r_pp - perp2[intersects])
            )
            d2 = np.maximum(0.0, d2 - chord)

    # 4) Threshold on the *corrected* distance: keep pairs whose radial
    #    separation (after size correction) lies within rmax.
    keep = d2 * d2 <= r2max
    if keep.any():
        iout = (i_idx[keep] + 1).tolist()
        jout = (j_idx[keep] + 1).tolist()
        dout = d2[keep].tolist()
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
