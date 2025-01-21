import numpy as np
import math
import pywt
import logging, sys
import scipy.stats as stats
import scipy.optimize

def hnorm(x, weights = None):
    # x is a numpy array
    num_rows, num_cols = x.shape
    if not weights:
        weights = [1]*num_rows
    n = sum(weights)
    sd = []
    for col in x.T:
        sd.append(math.sqrt(pywt.dwt(col)))
        # sd <- sqrt(apply(x, 2, wvar, w = weights))
    if num_cols == 2:
        hh = sd*(1/n)^(1/6)
    elif num_cols == 3:
        hh = sd * (4/(5 * n))^(1/7)
    else:
        logging.error("data with >3 dimensions are not allowed.")
        sys.exit(2)
    return hh

def phi6(x):
    return (np.power(x,6) - 15 * np.power(x,4) + 45 * np.power(x,2) - 15) * stats.norm.pdf(x)


def phi4(x):
    return (np.power(x,4) - 6 * np.power(x,2) + 3) * stats.norm.pdf(x)

def sj(x, h):
    # Tested function
    n = len(x)
    lambda_ = np.percentile(x, 75) - np.percentile(x, 25)
    a = 0.92 * lambda_ * n^(-1/7)
    b = 0.912 * lambda_ * n^(-1/9)
    W = W = np.reshape([element for element in x for i in range(n)],(n,n))
    W = W - np.full((n,n),x,order='C')
    W1 = phi6(W/b)
    tdb = np.matmul(np.matmul([1]*n,W1),[1]*n)
    tdb = -tdb / (n * (n - 1) * math.pow(b,7))

    W1 = phi4(W/a)
    sda = np.matmul(np.matmul([1]*n,W1),[1]*n)
    sda = sda / (n * (n - 1) * math.pow(a,5))

    alpha2 = 1.357 * math.pow((abs(sda/tdb)),(1/7)) * math.pow(h,(5/7))
    W1 = phi4(W/alpha2)
    sdalpha2 = np.matmul(np.matmul([1]*n,W1),[1]*n)
    sdalpha2 = sdalpha2 / (n * (n - 1) * math.pow(alpha2,5))

    result = math.pow((stats.norm.pdf(0,scale=math.sqrt(2)) / (n * abs(sdalpha2))),0.2) - h
    return result

def hsj(x):
    h0 = hnorm(x)
    v0 = sj(x, h0)
    if v0>0:
        hstep = 1.1
    else:
        hstep = 0.9
    h1 = h0 * hstep
    v1 = sj(x, h1)
    while (v1 * v0 > 0):
        h0 = h1
        v0 = v1
        h1 = h0 * hstep
        v1 = sj(x, h1)
    return h0 + (h1 - h0) * abs(v0)/(abs(v0) + abs(v1))
  

def h_select(x, y = None, weights = None, group = None, nbins=None, ndim=0, density=None, nobs=0):
#   source : https://github.com/cran/sm/blob/master/R/hselect.r
#   group   <- as.numeric(factor(group))
#   data    <- sm.check.data(x, y, weights = weights, group = group, ...)
#   x       <- data$x
#   y       <- data$y
#   weights <- data$weights
#   group   <- data$group
#   nobs    <- data$nobs
#   ndim    <- data$ndim
#   density <- data$density
#   opt     <- data$options
    if group:
        h_all = []
        for igroup in range(0, len(group),1):
            if ndim == 1:
                h_igroup = h_select(x[igroup], y[igroup],weights[igroup])
            else:
                h_igroup = h_select(x[igroup,], y[igroup,],weights[igroup,])
            h_all.append(h_igroup)
        h_all_len = len(h_all)
        h_all = np.array(h_all)
        h_all.reshape((h_all_len/ndim,ndim))
        h_gmean = []
        for col in h_all.T:
            sum = 0
            for item in col:
                sum+= math.log(item)
            h_gmean.append(math.exp(sum/len(col)))
        return h_gmean
    
    if ndim == 1:
        df = 6
    elif ndim == 1:
        df = 12

    if density:
        method = "normal"
    else:
        method = "df"

    if ndim==3 and not (density and method=="normal"):
        logging.error("bandwidth selection not available for 3 dimensions.")
        sys.exit(2)

    if not nbins:
        nbins=round((nobs > 100) * 8 * math.log(nobs) / ndim)

    data = [{"x", x}, {"means", y}, {"x.freq", [1]*nobs}, {"devs", [0]*nobs}, {"hweights",[1]*nobs}]
    sd = math.sqrt(np.diag(x))
    start = sd / 2
    data.update[{"sd", sd}]

    if (density & method == "normal"):
        return(hnorm(x))
    elif (density & method == "sj"):
        if (ndim > 1):
            logging.error("Sheather-Jones method requires 1-d data")
            sys.exit(2)
        return hsj(x)
    else:
        if (density):
            crit_type = "dens"
        else:
            crit_type = "reg"
        fname = method+".crit."+crit_type
        result = scipy.optimize.minimize(fname, [math.log(start/8), math.log(start*4)])
        h_result = math.exp(min(result))
    if ndim == 2:
        h_result = h_result*2
    return h_result


def sm_density_1d(x, h, eval_points, model = "none", hmult=1, weights = None, group = None):
    # est <- sm.density.eval.1d(x, h, weights)
    if "weights" in h.columns:
        h["weights"] = [1]*len(x)
    xlim = [min(x)-(range(x))/4, max(x)+(range(x))/4]
    ngrid=100
    if not eval_points:
        eval_points = list(range(int(xlim[0]),int(xlim[1])+1,1))
    xnew = eval_points
    n= len(x)
    neval = len(xnew)
    rep = [ele for ele in xnew for i in range(n)]
    W = np.array(rep)
    W.reshape((neval,n))
    hweights = h["weights"]
    W1 = np.array([ele for ele in hweights for i in range(neval)])
    W1.reshape((neval,n))
    W = math.exp(-0.5*(W/(hmult*h*W1))^2)/W1
    est = np.matmul(W,weights/(sum(weights) * math.sqrt(2 * math.pi) * hmult * h))
    # return [xnew, est, h*hmult, hweights, weights]
    return est


#     replace.na(opt, h.weights, rep(1, length(x)))
#     replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) +
#         diff(range(x))/4))
#     replace.na(opt, ngrid, 100)
#     hmult <- opt$hmult
#     h.weights <- opt$h.weights
#     xlim <- opt$xlim
#     ngrid <- opt$ngrid
#     replace.na(opt, eval.points, seq(xlim[1], xlim[2], length = ngrid))
#     xnew <- opt$eval.points
#     n <- length(x)
#     neval <- length(xnew)
#     W <- matrix(rep(xnew, rep(n, neval)), ncol = n, byrow = TRUE)
#     W <- W - matrix(rep(x, neval), ncol = n, byrow = TRUE)
#     W1 <- matrix(rep(h.weights, neval), ncol = n, byrow = TRUE)
#     W <- exp(-0.5 * (W/(hmult * h * W1))^2)/W1
#     est <- W %*% weights/(sum(weights) * sqrt(2 * pi) * hmult * h)
#     invisible(list(eval.points = xnew, estimate = as.vector(est),
#         h = h * hmult, h.weights = h.weights, weights = weights))

def sm_density_2d(x, y, h, xnew=None, ynew=None, hmult=1, eval_type = "points", weights = None):
    xlim = range(x)
    ylim = range(y)
    ngrid = 50
    if "weights" in h.columns:
        h["weights"] = [1]*len(x)
    if not xnew:
        xnew = list(range(int(xlim[0]),int(xlim[1])+1,(xlim[1]-xlim[0])/ngrid))
    if not ynew:
        ynew = list(range(int(ylim[0]),int(ylim[1])+1,(ylim[1]-ylim[0])/ngrid))
    n = len(x)
    nnew = len(xnew)
    hweights = h["weights"]
    W1 = np.array([ele for ele in xnew for i in range(n)])
    W1.reshape((nnew,n))
    W1_minus = np.array([x]*nnew)
    W1_minus.reshape((nnew//n,n))
    W1 = W1-W1_minus
    W2 = np.array([ele for ele in hweights for i in range(nnew)])
    W2.reshape((nnew,n))
    Wx = math.exp(-0.5 * (W1/(hmult * h[1] * W2))^2)/W2
    Wy = math.exp(-0.5 * (W1/(hmult * h[2] * W2))^2)/W2
    if (eval.type == "points"):
        est = np.matmul((Wx * Wy), weights)/(sum(weights) *
            2 * math.pi * h[1] * h[2] * hmult^2)
    else:
        est = np.matmul(Wx, (weights * np.transpose(Wy)))/(sum(weights) * 2 *
            math.pi * h[1] * h[2] * hmult^2)
    return est
    #  opt <- sm.options(options)
    # replace.na(opt, xlim, range(x))
    # replace.na(opt, ylim, range(y))
    # replace.na(opt, ngrid, 50)
    # replace.na(opt, h.weights, rep(1, length(x)))
    # if (missing(xnew))
    #     xnew <- seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid)
    # if (missing(ynew))
    #     ynew <- seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid)
    # n <- length(x)
    # nnew <- length(xnew)
    # h.weights <- opt$h.weights
    # hmult <- opt$hmult
    # W1 <- matrix(rep(xnew, rep(n, nnew)), ncol = n, byrow = TRUE)
    # W1 <- W1 - matrix(rep(x, nnew), ncol = n, byrow = TRUE)
    # W2 <- matrix(rep(h.weights, nnew), ncol = n, byrow = TRUE)
    # Wx <- exp(-0.5 * (W1/(hmult * h[1] * W2))^2)/W2
    # W1 <- matrix(rep(ynew, rep(n, nnew)), ncol = n, byrow = TRUE)
    # W1 <- W1 - matrix(rep(y, nnew), ncol = n, byrow = TRUE)
    # Wy <- exp(-0.5 * (W1/(opt$hmult * h[2] * W2))^2)/W2
    # if (eval.type == "points")
    #     est <- as.vector(((Wx * Wy) %*% weights)/(sum(weights) *
    #         2 * pi * h[1] * h[2] * hmult^2))
    # else est <- (Wx %*% (weights * t(Wy)))/(sum(weights) * 2 *
    #     pi * h[1] * h[2] * hmult^2)
    # invisible(list(eval.points = cbind(xnew, ynew), estimate = est,
    #     h = h * hmult, h.weights = h.weights, weights = weights))


def sm_density(X, eval_points, model = "none", weights = None, group = None, display="none", nobs=0, ndim=0):
    # source : https://github.com/cran/sm/blob/master/R/density.r
    nbins = round((nobs > 500) * 8 * math.log(nobs) / ndim)
    rawdata = [{"nbins", nbins}, {"x", x}, {"nobs", nobs}, {"ndim", ndim}]
    h = h_select(x, y = None, weights = weights, nbins=0)
    if X.ndim == 1:
        # replace.na(opt, xlab, x.name)
        # replace.na(opt, ylab, "Probability density function")
        est = sm_density_1d(x, h, eval_points, model, weights, rawdata)

    elif X.ndim == 2:
        x = X[0]
        y = X[1]
        est = sm_density_2d(x, y, h, eval_points, model, weights, rawdata)
    return est
