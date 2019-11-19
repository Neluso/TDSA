from numpy.fft import rfft, rfftfreq
from numpy import argmax, where, ones, abs  # array functions
from numpy import polyfit


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


def deconv():
    return 0
