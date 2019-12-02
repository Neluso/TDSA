from numpy import *
from scipy.optimize import minimize
from TDS_constants import *


def transfer_function(n, k, L, f):
    n_quo = 4 * n * n_aire / (n + n_aire)**2
    exp_term = exp(-1j * (n - n_aire) * (2 * pi * L) / c_0)
    fp = 1 / (1 - ())
    return 0