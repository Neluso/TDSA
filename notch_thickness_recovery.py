from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
from matplotlib.pyplot import *
from scipy.stats import cauchy
from scipy.signal import windows, filter_design


def thickness_val(n, order_m, f_ref, f_sam):
    return c_0 * (order_m + 0.5) / (n - n_air) * (1 / f_ref - 1 / f_sam)


def thickness_sigma2(n, sig_n, order_m, f_ref, sig_f_ref, f_sam, sig_f_sam):  # sigs in relative error (0 to 1)
    aux1 = ((1 / f_ref - 1 / f_sam) * c_0 * (order_m + 0.5) / (n - n_air) ** 2) ** 2 * (sig_n * n) ** 2
    aux2 = (c_0 * (order_m + 0.5) / (n - n_air)) ** 2 * ((sig_f_ref / f_ref) ** 2 + (sig_f_sam / f_sam) ** 2)
    return aux1 + aux2


def thickness_sigma2_aprox(n, order_m, f_ref, sig_f_ref, f_sam, sig_f_sam):  # sigs in relative error (0 to 1)
    aux = (c_0 * (order_m + 0.5) / (n - n_air)) ** 2 * ((sig_f_ref / f_ref) ** 2 + (sig_f_sam / f_sam) ** 2)
    return aux


def objetive_func(x, loc, scale, a, m):
    return cauchy.pdf(x, loc, scale) + a * x + m


def cauchy_fit(params, *args):
    loc, scale, a, m = params
    x, y = args
    return sum((cauchy.pdf(x, loc, scale) + a * x + m - y)**2)


# t_ref, E_ref = read_slow_data('./data/notch/ref_notch.txt')
# t_sam, E_sam = read_slow_data('./data/notch/sam_notch.txt')
# t_ref, E_ref = read_1file('./data/marca_autodestructiva/ref.txt')
# t_ref, E_ref = read_1file('./data/marca_autodestructiva/edge.txt')
# t_sam, E_sam = read_1file('./data/marca_autodestructiva/sam.txt')
# t_sam, E_sam = read_1file('./data/marca_autodestructiva/edge.txt')
t_ref, E_ref = read_slow_data('./data/notch/barrido 200 a 300 paso 0.1.txt')
t_sam, E_sam = read_slow_data('./data/notch/barrido 300a 200 paso -0.1.txt')


# win_idx = int(round(E_ref.size / 2))
# window = windows.tukey(win_idx, 0.2)
# window = zero_padding(window, 0, E_ref.size - win_idx)
# window2 = windows.tukey(win_idx)
# window2 = zero_padding(window2, 0, E_ref.size - win_idx)
# E_ref *= window
# E_sam *= window


# delta_t_ref = mean(diff(t_ref))
# enlargement = 0 * E_ref.size
# E_ref = zero_padding(E_ref, 0, enlargement)
# t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
# E_sam = zero_padding(E_sam, 0, enlargement)
# t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
#
# f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, E_ref.size)
# f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, E_sam.size)
# # H_w = E_sam_w / E_ref_w
# w_filt = wiener_filter(E_ref_w, beta=noise_floor(f_ref, E_ref_w, 1))


figure(1)
plot(t_ref, E_ref, lw=1)
plot(t_sam, E_sam, lw=1)

# figure(2)
# plot(f_ref, toDb(E_ref_w), lw=1)
# plot(f_ref, toDb(E_ref_w * w_filt), lw=1)


show()
