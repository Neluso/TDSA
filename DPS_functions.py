from numpy.fft import rfft, rfftfreq
from numpy import argmax, where, ones, abs, zeros, concatenate  # array functions
from numpy import polyfit
from scipy import signal


def fourier_analysis(t_data, E_data, nSamp):
    samp_int = (t_data[1] - t_data[0])  # seconds
    E_data_w = rfft(E_data, n=nSamp)
    f_data = rfftfreq(nSamp, d=samp_int)  # Hz
    return f_data, E_data_w


def centre_loc(E_data):  # finds the centre of the pulse based on... TODO improve algorithm
    t_0_pos = argmax(abs(E_data))
    return t_0_pos


def noise_floor(freq, E_data, f_lim):
    f_lim_idx = where(freq >= f_lim)[0][0]
    p = polyfit(freq[f_lim_idx:], E_data[f_lim_idx:], 1)
    return p[1] * ones(freq.size)


def deconv():  # TODO make function
    return 0


def zero_padding(t_val, E_val, n0_l, n0_r):
    t_zpd = concatenate((zeros(n0_l), t_val, zeros(n0_r)))
    E_zpd = concatenate((zeros(n0_l), E_val, zeros(n0_r)))
    return t_zpd, E_zpd


def force_exp_windowing():  # TODO return
    return 0


def bh_windowing(t_val, E_val):
    win = signal.blackmanharris(E_val.size)
    center_val = centre_loc(E_val)
    center_win = centre_loc(win)
    delta_center = center_win - center_val
    t_win, E_win = zero_padding(t_val, E_val, delta_center, 0)
    t_win = t_win[:E_val.size]
    E_win = E_win[:E_val.size] * win
    return t_win, E_win


def hann_windowing():
    return 0
