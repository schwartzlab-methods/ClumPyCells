import logging, sys
import numpy as np
import math
from scipy.stats import norm
from scipy.interpolate import interp1d


def unnormdensity(x, weights=None, from_=None, to_=None, num=None):
    weights = np.array(weights)
    w = weights
    totmass = np.sum(w)
    out = density(x, weights=weights, from_=from_, to_=to_, n=num)

    return_val = np.array(out)
    return return_val


def bindist(x, weights, lo, up, n):
    # https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/src/massdist.c
    ixmin = 0
    ixmax = n - 2
    xdelta = (up - lo) / (n - 1)
    y = [0] * 2 * n
    weights = weights.squeeze()
    for i in range(0, len(x), 1):
        xpos = (x[i] - lo) / xdelta
        ix = int(np.floor(xpos))
        fx = xpos - ix
        wi = weights[i]
        if ixmin <= ix and ix <= ixmax:
            y[ix] += (1 - fx) * wi
            y[ix + 1] += fx * wi
        elif ix == -1:
            y[0] += fx * wi
        elif ix == ixmax + 1:
            y[ix] += (1 - fx) * wi
    return y


def bandwidth(x):
    x = np.array(x)
    if len(x) < 2:
        logging.error("At least 2 points are needed to calculate the bandwidth")
    hi = np.std(x)
    q75, q25 = np.percentile(x, [75, 25])
    iqr = q75 - q25
    if hi < iqr / 1.34:
        lo = hi
    else:
        lo = iqr / 1.34
    return 0.9 * lo * len(x) ** (-0.2)


def density(
    x, bw=None, adjust=1, kernel="gaussian", weights=None, n=512, from_=None, to_=None
):
    # https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/density.R

    if bw == None:
        bw = bandwidth(x)
    nx = len(x)
    if len(weights) != nx:
        logging.error("'x' and 'weights' have unequal length")
        sys.exit(2)
    if np.any(weights < 0):
        logging.error("'weights' must not be negative")
        sys.exit(2)
    wsum = np.sum(weights)
    if wsum != 1:
        logging.warning("sum(weights) != 1  -- will not get true density")

    n_user = n
    n = max(n, 512)
    if n > 512:
        n = math.pow(2, math.ceil(math.log2(n)))
    bw = adjust * bw
    if bw <= 0:
        logging.error("'bw' is not positive.")
        sys.exit(2)

    if from_ == None:
        from_ = min(x) - 3 * bw
    if to_ == None:
        to_ = max(x) + 3 * bw

    lo = from_ - 4 * bw
    up = to_ + 4 * bw

    step = (up - lo) / n
    kords = list(np.arange(0, 2 * (up - lo), step))
    n = int(n)
    for i in range(n + 1, 2 * n, 1):
        kords[i] = kords[2 * n - i]

    if kernel == "gaussian":
        kords = norm.pdf(kords, scale=bw)
    elif kernel == "epanechnikov":
        a = bw * math.sqrt(5)
        ax = abs(kords)
        if ax < a:
            kords = 3 / 4 * (1 - (ax / a) ^ 2) / a
        else:
            kords = 0
    # This bins weighted distances
    y = np.array(bindist(x, weights, lo, up, n))

    fft_y = np.fft.fft(y)
    fft_kords = np.fft.fft(kords)
    z = fft_y * np.conjugate(fft_kords)
    kords = np.fft.ifft(z)
    kords = np.maximum(0, np.real(kords))[0:n]
    xords = list(np.arange(lo, up, step))

    y = interp1d(xords, kords, fill_value="extrapolate")

    step = (to_ - from_) / n_user
    x = list(np.arange(from_, to_, step))
    result = np.real(y(x))
    return result
