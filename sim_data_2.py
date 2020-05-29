from TDSA import *
from numpy.fft import *
from numpy.fft import fft as fft_func
from time import time_ns
import lmfit


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layer_number = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
n_subs = 1.25  # - 0.000 * 1j  # substrate refractive index -- cork 2.0


class Layer:
    def __init__(self, refactive_index, extintion_coefficient, layer_thickness, angle, freq):
        self.n = refactive_index  # may be an array depending from frequency
        self.k = extintion_coefficient  # must be an array depending on frequency if refractive index is it too
        self.complex_n = refactive_index - 1j * extintion_coefficient
        self.d = layer_thickness
        self.theta = angle
        self.phase = - 1j * 2 * pi * self.complex_n * (2 * self.d * cos(self.theta)) * freq / c_0


def H_sim(freq, layer_defs):
    r_coeffs = ones((layer_number+2, layer_number+2, freq.size), dtype=complex)
    t_coeffs = ones((layer_number+2, layer_number+2, freq.size), dtype=complex)
    for i in range(layer_number+2):
        for j in range(layer_number+2):
            n_cos_i = layer_defs[i].complex_n * cos(layer_defs[i].theta)
            n_cos_j = layer_defs[j].complex_n * cos(layer_defs[j].theta)
            denom_ij = (n_cos_i + n_cos_j)
            t_coeffs[i, j] = 2 * n_cos_i / denom_ij
            r_coeffs[i, j] = (n_cos_i - n_cos_j) / denom_ij
    print(real(t_coeffs[:, :, 0]))
    print(real(r_coeffs[:, :, 0]))
    # quit()
    H_i = r_coeffs[layer_number, layer_number+1]
    for i in -sort(-arange(1, layer_number+1)):
        exp_beta = exp(layer_defs[i].phase)
        H_i = r_coeffs[i-1, i] + (t_coeffs[i-1, i]*t_coeffs[i, i-1] * H_i * exp_beta) / (1 + r_coeffs[i-1, i] * H_i * exp_beta)
    return exp(layer_defs[0].phase) * H_i


# t_ref, E_ref = read_1file('./data/sim_resources/transmision_ref.txt')
t_ref, E_ref = read_1file('./data/sim_resources/refletion_ref.txt')
t_ref2, E_ref2 = read_1file('./data/sim_resources/ref_reflexion.txt')


t_ref *= 1e-12
t_ref2 *= 1e-12
# print(t_ref)
# quit()

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# f_ref2, E_ref_w2 = fourier_analysis(t_ref2, E_ref2)
# f_min, f_max = f_min_max_idx(f_ref, 0.2, 1.1)
layers = list()

f_ref *= 1e-12
dair = 0  # 20e-8
medium = Layer(n_air, 0, dair, deg_in * pi / 180, f_ref*1e12)
layers.append(medium)

n1 = 3 * ones(f_ref.size)
k1 = 1 * f_ref
d1 = 50e-6
layer_1 = Layer(n1, k1, d1, arcsin(snell_sin / n1), f_ref*1e12)
layers.append(layer_1)

n2 = 1.5 * ones(f_ref.size)
k2 = 1 * f_ref
d2 = 50e-6
layer_2 = Layer(n2, k2, d2, arcsin(snell_sin / n2), f_ref*1e12)
layers.append(layer_2)

n3 = 3 * ones(f_ref.size)
k3 = 1 * f_ref
d3 = 50e-6
layer_3 = Layer(n3, k3, d3, arcsin(snell_sin / n3), f_ref*1e12)
layers.append(layer_3)

substrate = Layer(n_subs, 0.0, 0, arcsin(snell_sin / n_subs), f_ref*1e12)
layers.append(substrate)

f_ref *= 1e12


H_0 = H_sim(f_ref, layers)
wh = open('H_0.txt','w')
for i in range(f_ref.size):
    wh.write(str(f_ref[i])+','+str(real(H_0[i]))+','+str(imag(H_0[i]))+'\n')
# quit()
figure(20)
plot(f_ref, abs(H_0))
xlim([0, 1e12])
figure(21)
plot(f_ref, unwrap(angle(H_0)))
xlim([0, 1e12])
# E_sim = irfft(rfft(E_ref) * H_0)
E_sim = ifft(E_ref_w * H_0)

figure(1)
plot(t_ref, E_ref, lw=1, label='ref1')
plot(t_ref, E_sim, lw=1, label='sim1')
legend()
xlim([t_ref[0], t_ref[-1]])


# layers = list()
#
# f_ref2 *= 1e-12
# dair = 0  # 20e-9
# medium = Layer(n_air, 0, dair, deg_in * pi / 180, f_ref2*1e12)
# layers.append(medium)
#
# n1 = 2.5 * ones(f_ref2.size)
# k1 = 0.01 * f_ref2
# d1 = 50e-6
# layer_1 = Layer(n1, k1, d1, arcsin(snell_sin / n1), f_ref2*1e12)
# layers.append(layer_1)
#
# n2 = 2.8 * ones(f_ref2.size)
# k2 = 0.01 * f_ref2
# d2 = 50e-6
# layer_2 = Layer(n2, k2, d2, arcsin(snell_sin / n2), f_ref2*1e12)
# layers.append(layer_2)
#
# n3 = 2.7 * ones(f_ref2.size)
# k3 = 0.01 * f_ref2
# d3 = 50e-6
# layer_3 = Layer(n3, k3, d3, arcsin(snell_sin / n3), f_ref2*1e12)
# layers.append(layer_3)
#
# substrate = Layer(n_subs, 0.0, 0, arcsin(snell_sin / n_subs), f_ref2*1e12)
# layers.append(substrate)
#
# f_ref2 *= 1e12
#
#
# H_02 = H_sim(f_ref2, layers)
# figure(20)
# plot(f_ref2, abs(H_02))
# xlim([0, 2e12])
# figure(21)
# plot(f_ref2, unwrap(angle(H_02)))
# xlim([0, 2e12])
# E_sim2 = irfft(rfft(E_ref2) * H_02)
# # E_sim2 = ifft(E_ref_w2 * H_02)
#
#
# figure(1)
# plot(t_ref2, E_ref2, lw=1, label='ref2')
# plot(t_ref2, E_sim2, lw=1, label='sim2')
# legend()
# xlim([t_ref2[0], t_ref2[-1]])

show() 
