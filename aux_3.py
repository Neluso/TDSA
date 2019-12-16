from read_data import *
from scipy.signal import windows
from DSP_functions import *
from jepsen_index import jepsen_index
from aux_functions import *


t_ref, E_ref = read_slow_data('./data/cork/ref_cork.txt')
t_sam, E_sam = read_slow_data('./data/cork/sam_cork.txt')

# t_ref, E_ref = read_1file('./data/demo_data/test_ref.txt')
# t_sam, E_sam = read_1file('./data/demo_data/test_sam.txt')
print(t_ref)

t_ref -= t_ref[0]
t_sam -= t_sam[0]

t_ref *= 1e-12
t_sam *= 1e-12
thick = 2.6e-3
print(t_ref)



wind_idx = int(round(E_ref.size / 2))
zero_idx = E_ref.size - wind_idx
window = windows.tukey(wind_idx, 0.3)
window = zero_padding(window, int(round(zero_idx / 3)), int(round(zero_idx * 2 / 3)))
E_ref *= window
E_sam *= window

nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)
n_f, alpha_f, n_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, thick)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
H_w = E_sam_w / E_ref_w


t_ref *= 1e12
t_sam *= 1e12
f_ref *= 1e-12
f_sam *= 1e-12

figure(10)
plot(f_ref, abs(H_w), lw=1)
xlim([0, 2])
figure(11)
plot(f_ref, unwrap(angle(H_w)), lw=1)
xlim([0, 2])

figure(1)
plot(t_ref, E_ref, lw=1)
plot(t_sam, E_sam, lw=1)
plot(t_ref, window, lw=1)

figure(2)
plot(f_ref, toDb(E_ref_w), lw=1)
plot(f_sam, toDb(E_sam_w), lw=1)
xlim([0, 2])

figure(3)
title('Refractive index')
plot(f_ref, n_f, lw=1)
xlim([0, 2])

figure(4)
title('Absorption coeff')
plot(f_ref, 0.01 * alpha_f, lw=1)
xlim([0, 2])


show()
