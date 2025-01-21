import logging
import math
import os
import pickle
import sys
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .breakpts import *
from .closepairs import *
from .pointPattern import *
from .sm_density import *
from .unnormdensity import *
from .window import *


def convertToWindow(domain):
    # assume domain[0] is a list for x, domain[1] is a list for y
    return window(domain[0], domain[1])


# TO DO: convert all lists into numpy
def pmin(a, b):
    return np.minimum(a, b)


def hang(d, r):
    nr, nc = r.shape
    d = np.array(d).reshape((nr, nc))
    answer = np.zeros((nr, nc))
    mask = d < r
    answer[mask] = np.arccos(d[mask] / r[mask])
    return answer


def edgecorrection(points, r, w=None, maxweight=100):
    if w is None:
        w = points.window
    x = points.getX()
    y = points.getY()
    n = len(x)
    nrow, ncol = r.shape
    assert nrow == n, "number of row should match the number of points"
    if nrow * ncol == 0:
        return r

    xrange = w.getXrange()
    yrange = w.getYrange()
    dL = x - xrange[0]
    dR = xrange[1] - x
    dD = y - yrange[0]
    dU = yrange[1] - y

    dL_reshape = dL.reshape(-1, 1)
    dR_reshape = dR.reshape(-1, 1)
    dD_reshape = dD.reshape(-1, 1)
    dU_reshape = dU.reshape(-1, 1)

    bLU = np.arctan2(dU_reshape, dL_reshape)
    bLD = np.arctan2(dD_reshape, dL_reshape)
    bRU = np.arctan2(dU_reshape, dR_reshape)
    bRD = np.arctan2(dD_reshape, dR_reshape)
    bUL = np.arctan2(dL_reshape, dU_reshape)
    bUR = np.arctan2(dR_reshape, dU_reshape)
    bDL = np.arctan2(dL_reshape, dD_reshape)
    bDR = np.arctan2(dR_reshape, dD_reshape)

    aL = hang(dL_reshape, r)
    aR = hang(dR_reshape, r)
    aD = hang(dD_reshape, r)
    aU = hang(dU_reshape, r)

    cL = pmin(aL, bLU) + pmin(aL, bLD)
    cR = pmin(aR, bRU) + pmin(aR, bRD)
    cU = pmin(aU, bUL) + pmin(aU, bUR)
    cD = pmin(aD, bDL) + pmin(aD, bDR)

    ext = cL + cR + cU + cD
    ext = ext.reshape(nrow, ncol)

    weight = 1 / (1 - ext / (2 * math.pi))
    weight = weight.squeeze()
    return weight


def edgetrans(X, Y=None, W=None, exact=False, paired=False, dx=None, dy=None):
    # https://github.com/spatstat/spatstat.core/blob/a560307cd0a8e4f083ed4138d2bb968314da02e0/R/edgeTrans.R
    if not Y:
        Y = X
    if not W:
        W = X.getWindow()
    if not dx or not dy:
        xX = X.getX()
        xY = Y.getX()
        nX = len(xX)
        nY = len(xY)
        if paired:
            if nX != nY:
                logging.error("X and Y should have equal length when paired=True")
                sys.exit(2)
            dx = np.subtract(np.array(xX), np.array(xY))
            dy = np.subtract(np.array(X.getY()), np.array(Y.getY()))
        else:
            dx = np.subtract.outer(np.array(xX), np.array(xY))
            dy = np.subtract.outer(np.array(X.getY()), np.array(Y.getY()))
    else:
        # dx, dy given
        if paired:
            if dx != dy:
                logging.error("len(dx) != len(dy)")
                sys.exit(2)
        else:
            # dx, dy are matrices
            nX = len(dx)
            nY = len(dx[0])
    # For irregular polygons, exact evaluation is very slow;
    # so use pixel approximation, unless exact=True
    if W.getType() == "rectangular":
        # Fast code for this case
        wide = np.diff(W.getXrange())
        high = np.diff(W.getYrange())
        weight = wide * high / ((wide - abs(dx)) * (high - abs(dy)))
    elif W.getType() == "polygonal":
        n = len(dx)
        weight = []
        if n > 0:
            for i in range(0, n, 1):
                Wshift = window(dx[i], dy[i])
                weight.append(window.overlap(W, Wshift))
    if not paired:
        weight = np.array(weight).reshape((nX, nY))

    return weight


def implemented_for_K(correction, windowtype, explicit):
    if "best" in correction:
        # select best available correction
        if windowtype != "mask":
            for i in range(0, len(correction)):
                if correction[i] == "best":
                    correction[i] = "isotropic"
        else:
            for i in range(0, len(correction)):
                if correction[i] == "best":
                    correction[i] = "translate"
    else:
        if windowtype == "mask":
            iso = []
            for i in range(0, len(correction)):
                if correction[i] == "isotropic":
                    iso.append(1)
                    logging.warning(
                        "Isotropic correction not implemented for binary masks"
                    )
                else:
                    iso.append(0)
            if sum(iso) == len(correction) and explicit:
                logging.error(
                    "Isotropic correction not implemented for binary masks, but got only isotropic as choice"
                )
                sys.exit(2)
            correction.remove("isotropic")
    return correction


## Smooth Estimate of Weighted Second Moment Density
def sewsmod(d, ff, wt, Ef, rvals, method="smrep", nwtsteps=500):
    # d is a list of interpoint distance
    # ff = f(M1, M2) where M1, M2 are marks at the two points
    # wt = edge correction weight

    if method == "smrep":
        fw = [ff[i] * wt[i] for i in range(len(ff))]
        nw = round(nwtsteps * wt / max(wt))
        drep_w = [distance * nw for distance in d]
        nfw = round(nwtsteps * fw / max(fw))
        drep_fw = [distance * nfw for distance in d]
        est = sm_density(drep_fw, eval_points=rvals, display="none", nbins=0)
        numerator = est * sum(fw) / sum(est)
        est0 = sm_density(drep_w, eval_points=rvals, display="none", nbins=0)
        denominator = est0 * (sum(wt) / sum(est0)) * Ef
        result = numerator / denominator
    elif method == "sm":
        fw = [ff[i] * wt[i] for i in range(len(ff))]
        est = sm_density(d, weights=fw, eval_points=rvals, display="none", nbins=0)
        numerator = est * sum(fw) / sum(est)
        est0 = sm_density(d, weights=wt, eval_points=rvals, display="none", nbins=0)
        denominator = est0 * (sum(wt) / sum(est0)) * Ef
        result = numerator / denominator
    elif method == "density":
        # ff = np.array(ff)np.multiply(ff,wt)#
        w = [ff[i] * wt[i] for i in range(len(ff))]

        Kf = unnormdensity(
            d, weights=w, from_=min(rvals), to_=max(rvals), num=len(rvals)
        )
        # smooth estimate of kappa_1

        K1 = unnormdensity(
            d, weights=wt, from_=min(rvals), to_=max(rvals), num=len(rvals)
        )
        result = Kf / (Ef * K1)
    else:
        print("currently not support loess")
    return result


def rmax_rule(W, area):
    # reference: Kest(X) at https://spatstat.org/FAQ.html
    ripley_rmax = 0.25 * W.minEdge()
    large_rmax = math.sqrt(1000 / (math.pi * area))
    return min(ripley_rmax, large_rmax)


def pickoption(
    key,
    keymap,
    what="option",
    exact=False,
    list_on_err=True,
    die=True,
    multi=False,
    allow_all=True,
):
    if len(key) == 0:
        logging.error("Argument " + key + " has length zero")
        sys.exit(2)
    key = set(key)
    if not multi and len(key) > 1:
        logging.error("Must specify only one " + what + " " + key)
        sys.exit(2)
    allow_all == allow_all and multi
    if allow_all and ("all" in key):
        id = keymap
    elif exact:
        id = []
        for k in key:
            if k in keymap:
                id.append(keymap[k])
    else:
        id = []
        pkeymap = str(keymap)
        for k in key:
            if k in pkeymap:
                id.append(keymap[k])
    key = []
    for i in id:
        key.append(keymap[i])
    return key


def visualize_result(funs):
    names, kmm = [], []
    for i in funs:
        names.append(i)
        kmm.append(funs[i][0])
        kmm.append(funs[i][1])
    kmm = np.asarray(kmm).T
    iterables = [names, ["iso", "trans"]]
    index = pd.MultiIndex.from_product(iterables)
    kmm = pd.DataFrame(kmm, columns=index)
    print(kmm)


def markcorr(
    X,
    r=None,
    correction=["translate", "isotropic"],
    method=["density"],
    normalise=True,
    weights=None,
    multitype=2,
    enable_type=False,
    saveImage=True,
    savefolder="./result/",
    remove_zeros=True,
    saveCache=True,
    pp=None,
):
    # obtain point pattern values
    markx = X.getMarks()
    npts = len(X.getX())
    W = X.getWindow()
    d = X.getD()
    if not weights:
        weights = np.ones(npts)
    else:
        if np.any(weights == 0):
            sys.exit(2)

    # determine the maximum range of r
    rmaxdefault = rmax_rule(W, npts / X.getWindow().getArea())
    breaks = handle_r_b(r, None, W, rmaxdefault=rmaxdefault)
    r = breaks.getR()
    rmax = breaks.getMax()

    r_out = pd.DataFrame({"r": r})
    r_out.to_csv(savefolder + "r.csv")

    if len(method) > 1:
        logging.error("Select only one method, please")
        sys.exit(2)
    if (method == "density") and ("even" in breaks.columns):
        logging.error("Evenly spaced r values are required if method=density")
        sys.exit(2)

    # available selection of edge corrections depends on window
    correction = pickoption(
        correction,
        {
            "none": "none",
            "border": "border",
            "bord.modif": "bord.modif",
            "isotropic": "isotropic",
            "Ripley": "isotropic",
            "translate": "translate",
            "translation": "translate",
            "best": "best",
        },
        what="correction",
        multi=True,
    )
    correction = implemented_for_K(correction, W.getType(), False)

    # TODO: enable type is meant for grouping rows based on the "type" array.
    # For now, the plan is to use a separate attribute in pointpattern
    if enable_type:
        dumi = pd.get_dummies(markx["type"], prefix=str(column) + "_")
        # for m in markx:
        #     if m
    for column in markx:
        if markx[column].dtype == "category":
            dumi = pd.get_dummies(markx[column], prefix=str(column))
            markx = markx.drop(columns=column)
            markx = markx.join(dumi)

            if remove_zeros:
                logging.warning(
                    "WARNING: Self correlation for catagorical values will be zero"
                )

    # enable multitype mark-cross correlation
    if multitype < 1:
        logging.error("number of types has to be greater than 0")
    if multitype > 2:
        colnames = markx.columns
        comb = combinations(colnames, multitype - 1)
        for i in list(comb):
            markx[str(i)] = markx[list(i)].mean(axis=1)

    # creating data structure to store the result
    funs = {}
    if remove_zeros == False:
        # find close pairs of points for all the points
        close = closepairs(X, rmax, d, pp=pp)
        dIJ = close["d"]
        I = close["i"]
        J = close["j"]
        XI = pointPattern(close["xi"], close["yi"], d, W)

    for coli in markx:
        for colj in markx:
            name = f"{coli} vs. {colj}"

            # Skip if the result is cached
            if os.path.exists(f"{savefolder}{name}.pkl"):
                continue

            result = []
            mari = markx[coli].to_numpy()
            marj = markx[colj].to_numpy()

            Ef = 0
            # Denominator
            Ef = np.mean(np.multiply(markx[coli], weights)) * np.mean(
                np.multiply(markx[colj], weights)
            )

            # check validity of denominator
            if Ef == 0:
                logging.error(
                    "Cannot normalise the mark correlation; the denominator is zero"
                )
                result.append([1] * len(r))
                result.append([1] * len(r))

            elif Ef < 0:
                logging.warning(
                    "Problem when normalising the mark correlation: the denominator is negative"
                )

            else:
                if remove_zeros:
                    # remove lines with all zeros
                    if coli != colj:
                        markpairs = markx[[coli, colj]]
                        validIndex = ~(markpairs[[coli, colj]] == 0).all(axis=1)
                        markpairs = markpairs[validIndex]
                        dpairs = d[validIndex]
                        mari = mari[validIndex]
                        marj = marj[validIndex]
                        xpairs = X.x[validIndex]
                        ypairs = X.y[validIndex]
                        X_temp = pointPattern(xpairs, ypairs, dpairs, W)

                    else:
                        validIndex = markx[coli] != 0
                        markpairs = markx.loc[validIndex]
                        dpairs = d[validIndex]
                        mari = mari[validIndex]
                        marj = marj[validIndex]
                        xpairs = X.x[validIndex]
                        ypairs = X.y[validIndex]
                        X_temp = pointPattern(xpairs, ypairs, dpairs, W)

                    if X_temp.n < 10 or np.mean(mari) == 0 or np.mean(marj) == 0:
                        result.append([1] * len(r))
                        result.append([1] * len(r))
                        name = coli + " vs. " + colj
                        funs[name] = result

                        if saveImage == True:
                            fig = plt.figure()
                            plt.plot(r, result[0], label="iso")
                            plt.plot(r, result[1], label="trans")
                            plt.plot(r, [1] * len(r), "--")
                            plt.xlabel("r")
                            plt.ylabel("kmm(r)")
                            plt.title(name)
                            plt.legend()
                            name = savefolder + coli + "_" + colj + ".png"
                            plt.savefig(fname=name)
                            plt.close("all")
                        continue

                    else:
                        # find close pairs of points
                        close = closepairs(X_temp, rmax, dpairs, pp=pp)
                        dIJ = close["d"]
                        I = close["i"]
                        J = close["j"]
                        XI = pointPattern(close["xi"], close["yi"], d, W)

                else:
                    nonZeroPoints = np.min(
                        [np.count_nonzero(markx[coli]), np.count_nonzero(markx[colj])]
                    )
                    if nonZeroPoints < 10:
                        if "isotropic" in correction:
                            result.append([1] * len(r))
                        if "translate" in correction:
                            result.append([1] * len(r))
                        name = coli + " vs. " + colj
                        funs[name] = result

                        if saveImage == True:
                            fig = plt.figure()
                            if len(correction) == 1:
                                if "isotropic" in correction:
                                    plt.plot(r, result[0], label="iso")
                                else:
                                    plt.plot(r, result[0], label="trans")
                            elif len(correction) == 2:
                                plt.plot(r, result[0], label="iso")
                                plt.plot(r, result[1], label="trans")
                            else:
                                plt.plot(r, result[0])
                            plt.plot(r, [1] * len(r), "--")
                            plt.xlabel("r")
                            plt.ylabel("kmm(r)")
                            plt.title(name)
                            plt.legend()
                            name = savefolder + coli + "_" + colj + ".png"
                            plt.savefig(fname=name)
                            plt.close("all")
                        continue
                # apply f to marks of close pairs of points
                mI = [mari[i - 1] for i in I]
                mJ = [marj[j - 1] for j in J]

                ff = np.multiply(mI, mJ)
                # ff = [mI[i] * mJ[i] for i in range(len(mI))]

                # Compute estimates
                if "none" in correction:
                    # uncorrected estimate
                    edgewt = [1] * len(dIJ)
                    # get smoothed estimate of mark covariance
                    Mnone = sewsmod(dIJ, ff, edgewt, Ef, r, method[0])
                    result.append(Mnone)

                if "translate" in correction:
                    XJ = pointPattern(close["xj"], close["yj"], W)
                    edgewt = edgetrans(XI, XJ, paired=True)
                    # get smoothed estimate of mark covariance
                    Mtrans = sewsmod(dIJ, ff, edgewt, Ef, r, method[0])
                    result.append(Mtrans)

                if "isotropic" in correction:
                    # Ripley isotropic correction
                    edgewt = edgecorrection(XI, np.array(dIJ).reshape(len(dIJ), 1))
                    # get smoothed estimate of mark covariance
                    Miso = sewsmod(dIJ, ff, edgewt, Ef, r, method[0])
                    result.append(Miso)

            if saveCache:
                f = open(f"{savefolder}{name}.pkl", "wb")
                pickle.dump(result, f)
                f.close()
            else:
                funs[name] = result

            if saveImage == True:
                fig = plt.figure()
                if len(correction) == 1:
                    if "isotropic" in correction:
                        plt.plot(r, result[0], label="iso")
                    else:
                        plt.plot(r, result[0], label="trans")
                elif len(correction) == 2:
                    plt.plot(r, result[0], label="iso")
                    plt.plot(r, result[1], label="trans")
                else:
                    plt.plot(r, result[0])
                plt.plot(r, [1] * len(r), "--")
                plt.xlabel("r")
                plt.ylabel("kmm(r)")
                plt.title(name)
                plt.legend()
                name = savefolder + coli + "_" + colj + ".png"
                plt.savefig(fname=name)
                plt.close("all")

    if saveCache:
        funs = {}
        for filename in os.listdir(savefolder):
            if filename.endswith(".pkl"):
                file_path = os.path.join(savefolder, filename)
                with open(file_path, "rb") as file:
                    try:
                        file_data = pickle.load(file)
                    except:
                        logging.ERROR(f"Missing File:{file_path}")
                        exit(1)
                funs[filename[:-4]] = file_data
    if saveImage:
        keys = list(funs.keys())
        ncol = len(markx.columns)
        fig, ax = plt.subplots(ncol, ncol, squeeze=False)

        for i in range(ncol):
            for j in range(ncol):
                key = keys[j + i * ncol]
                fun_curves = funs[key]

                # Plot the first curve (index 0) if it exists
                if len(fun_curves) > 0:
                    ax[i, j].plot(r, fun_curves[0], label="Curve 1", color="blue")

                # Plot the second curve (index 1) if it exists
                if len(fun_curves) > 1:
                    ax[i, j].plot(r, fun_curves[1], label="Curve 2", color="orange")

                # Plot the reference line
                ax[i, j].plot(r, [1] * len(r), "--", label="Reference", color="green")
                ax[i, j].set_title(key)

        for x in ax.flat:
            x.set(xlabel="r", ylabel="kmm")

        # Update the legend to handle one or two curves
        if len(correction) == 1:
            fig.legend([correction[0], "Reference"])
        else:
            fig.legend(["translate", "isotropic", "Reference"])

        fig.tight_layout()

        plt.savefig(savefolder + "all_images.png")
        plt.close("all")

    if saveCache:
        for filename in os.listdir(savefolder):
            if filename.endswith(".pkl"):
                file_path = os.path.join(savefolder, filename)
                try:
                    os.remove(file_path)
                except Exception as e:
                    logging.error(f"Error removing file: {file_path} - {e}")

    return r, funs
