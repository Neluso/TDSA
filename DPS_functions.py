from numpy.fft import rfft, rfftfreq
from numpy import exp, log, mean, diff  # numerical function
from numpy import argmax, where, abs, concatenate, flip  # array functions
from numpy import ones, zeros, linspace  # array creation functions
from matplotlib.pyplot import *
from numpy import polyfit
from scipy import signal


def fourier_analysis(t_data, E_data, nSamp):
    samp_int = mean(diff(t_data))  # seconds
    E_data_w = rfft(E_data, n=nSamp)
    f_data = rfftfreq(nSamp, d=samp_int)  # Hz
    return f_data, E_data_w


def centre_loc(E_data):  # finds the pulse centre  based on pulse maximum of absolute signal TODO improve algorithm
    t_0_pos = argmax(abs(E_data))
    return t_0_pos


def noise_floor(freq, E_data, f_lim):
    f_lim_idx = where(freq >= f_lim)[0][0]
    p = polyfit(freq[f_lim_idx:], E_data[f_lim_idx:], 1)
    return p[1] * ones(freq.size)


def zero_padding(val, n0_l, n0_r):
    if n0_l > 0:
        val = concatenate((zeros(n0_l), val))
    if n0_r > 0:
        val = concatenate((val, zeros(n0_r)))
    return val


def force_exp_func(data_size, top_half_size, decay_length, decay_minimum=1):  # decay_minimum default is 1%
    decay_rate = - log(decay_minimum / 100) / decay_length
    x = linspace(0, decay_length, num=decay_length)
    decay = exp(- x * decay_rate)
    half_point = int(data_size / 2)
    y = zeros(data_size)
    for i in range(y.size):
        if half_point - top_half_size <= i < half_point + top_half_size:
            y[i] = 1
        if half_point + top_half_size <= i < decay.size + half_point + top_half_size - 1:
            y[i] = decay[i - half_point - top_half_size]
    return y


def force_exp_windowing(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= force_exp_func(E_val.size, int(center_val_idx / 2), int(center_val_idx / 4))
    E_sub *= force_exp_func(E_sub.size, int(center_sub_idx / 2), int(center_sub_idx / 4))
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub


def bh_windowing(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= signal.blackmanharris(E_val.size)
    E_sub *= signal.blackmanharris(E_sub.size)
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub


def cheb_windowing(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= signal.chebwin(E_val.size, 80)
    E_sub *= signal.chebwin(E_sub.size, 80)
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub


def hann_windowing(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= signal.hann(E_val.size)
    E_sub *= signal.hann(E_sub.size)
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub


def tukey_windowing(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= signal.tukey(E_val.size)
    E_sub *= signal.tukey(E_sub.size)
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub
