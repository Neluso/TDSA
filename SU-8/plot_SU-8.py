from TDSA import *
from scipy.interpolate import Akima1DInterpolator


def nk_from_eps(e_w):
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


f_real, eps_real = read_1file('./real_fully_cross.csv')
f_imag, eps_imag = read_1file('./imag_fully_cross.csv')

e_r = Akima1DInterpolator(f_real, eps_real)
e_i = Akima1DInterpolator(f_imag, eps_imag)


n, k = nk_from_eps(e_r(f_real) + 1j * e_i(f_real))

plot(f_real, n)
show()
