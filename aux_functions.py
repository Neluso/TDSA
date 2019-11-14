from numpy import where  # array functions
from numpy import log10  # mathematical functions


def f_min_max_idx(freq):
    f_min_idx = 1
    f_max_idx = 100
    f_min = freq[f_min_idx]
    f_max = freq[f_max_idx]
    for frq in freq:
        if frq <= 1e11:  # 0.1 THz
            f_min = frq
        if frq <= 1e12:  # 1 THz
            f_max = frq
    f_min_idx = where(freq == f_min)[0][0]
    f_max_idx = where(freq == f_max)[0][0]
    return f_min_idx, f_max_idx


def toDb(x):
    return 20 * log10(abs(x))


def prettyfy(x, norm):
    return toDb(x / norm)


def nextpow2(i):
    n = 1
    while n < i:
        n *= 2
    return n


