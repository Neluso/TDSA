from numpy import array, sin, pi, concatenate
from numpy.fft import fft, fftfreq, ifft
import matplotlib.pyplot as plt

N = 2000
f = range(1, N + 1, 1)
f = 0.001 * array(f)
trans_f = sin(f + 1.5 * pi / 4)


for i in range(f.shape[0]):
    if trans_f[i] <= 0:
        trans_f[i] = 0


trans_t = ifft(trans_f)
bins = trans_t.shape[0]
mid_point = N - 100  # int(bins / 2)
front_trans = trans_t[:mid_point]
rear_trans = trans_t[mid_point:]
trans_t = concatenate((rear_trans, front_trans))


# plt.plot(f, trans_f)
t = 50 * f

exp_t = t[90:120]
exp_trans_t = trans_t[90:120]
exp_trans_f = fft(exp_trans_t)
exp_f = exp_t / 50  # fftfreq(exp_t.shape[0])
plt.plot(exp_f, exp_trans_f)
plt.show()