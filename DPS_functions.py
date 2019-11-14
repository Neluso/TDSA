from numpy.fft import rfft, rfftfreq
from numpy import argmax


def fourier_analysis(t_data, E_data, nSamp):
    samp_int = (t_data[1] - t_data[0])  # seconds
    E_data_w = rfft(E_data, n=nSamp)
    f_data = rfftfreq(nSamp, d=samp_int)  # Hz
    return f_data, E_data_w


def centre_loc(E_data):  # finds the centre of the pulse based on
    t_0_pos = argmax(abs(E_data))
    return t_0_pos


def noise_floor():
    return 0