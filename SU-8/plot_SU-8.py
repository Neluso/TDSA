from TDSA import *
from scipy.interpolate import Akima1DInterpolator, interp1d


def nk_from_eps(e_w):
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


f_real, eps_real = read_1file('./real_fully_cross.csv')
f_imag, eps_imag = read_1file('./imag_fully_cross.csv')
# plot(f_imag, eps_imag, '-')
# show()
# quit()

e_r = interp1d(f_real, eps_real, fill_value='extrapolate')
e_i = interp1d(f_imag, eps_imag, fill_value='extrapolate')

n, k = nk_from_eps(e_r(f_real) + 1j * e_i(f_real))
nh = open('./n_SU8.txt', 'w')
kh = open('./k_SU8.txt', 'w')

for i in range(f_real.size):
    nh.write(str(f_real[i]) + ',' + str(n[i]) + '\n')
    kh.write(str(f_real[i]) + ',' + str(k[i]) + '\n')
