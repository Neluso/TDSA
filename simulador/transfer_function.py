from simulador.sim_file import *


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
    H_i = r_coeffs[layer_number, layer_number+1]
    for i in -sort(-arange(1, layer_number+1)):
        exp_beta = exp(layer_defs[i].phase)
        H_i = r_coeffs[i-1, i] + (t_coeffs[i-1, i]*t_coeffs[i, i-1] * H_i * exp_beta) / (1 + r_coeffs[i-1, i] * H_i * exp_beta)
    return exp(layer_defs[0].phase) * H_i