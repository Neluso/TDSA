from TDS_constants import *
from DSP_functions import *
from read_data import read_1file
from aux_functions import *
from scipy.optimize import differential_evolution


def ct2(n_l, n_l_1):
    return 4 * n_l * n_l_1 / (n_l**2 + n_l_1**2)


def cr_l_1_l(n_l, theta_l, n_l_1, theta_l_1):  # from n_l-1 to n_l
    return (n_l_1 * cos(theta_l_1) - n_l * cos(theta_l)) / (n_l_1 * cos(theta_l_1) + n_l * cos(theta_l))


def theta_l(n_l, n_l_1, theta_l_1):
    return arcsin(n_l_1 * sin(theta_l_1) / n_l)


def phase_factor(n, k, thick, freq, theta_l):  # theta in radians
    n_comp = n + 1j * k
    omg = 2 * pi * freq
    thick *= cos(theta_l)
    return exp(2 * 1j * n_comp * omg * thick / c_0)


def H_sim(n_l, k_l, thick_l, freq):  # , theta):
    theta_l_1 = pi / 6  # theta
    n_l_1 = n_aire
    r_llp1 = 1
    theta_l_calc = theta_l(n_l, n_l_1, theta_l_1)
    exp_beta = phase_factor(n_l, k_l, thick_l, freq, theta_l_1)
    r_lm1l = cr_l_1_l(n_l, theta_l_calc, n_l_1, theta_l_1)
    tt = ct2(n_l, n_l_1)
    
    return r_lm1l + (tt * r_llp1 * exp_beta) / (1 + r_lm1l * r_llp1 * exp_beta)


# def objective_function(k, *args):
#     k0, k1, k2, k3, k4, k5, k6, k7 = k
#     E_val, E_meas, t_val = args
#     idx_k1 = centre_loc(E_val) - int(round(k1 / mean(diff(t_val))))
#     idx_k3 = centre_loc(E_val) - int(round(k3 / mean(diff(t_val))))
#     idx_k5 = centre_loc(E_val) - int(round(k5 / mean(diff(t_val))))
#     idx_k7 = centre_loc(E_val) - int(round(k7 / mean(diff(t_val))))
#     return k0 * roll(E_val, - idx_k1) + k2 * roll(E_val, - idx_k3) + k4 * roll(E_val, - idx_k5) + k6 * roll(E_val, - idx_k7)


def cost_function(k, *args):
    n1, k1, thick1 = k
    H_meas, freqs = args
    logs = log(abs(H_sim(n1, k1, thick1, freqs))) - log(abs(H_meas))
    angles = angle(H_sim(n1, k1, thick1, freqs)) - angle(H_meas)
    return sum(logs**2 + angles**2)


# Main script


t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/old/Paintmeter/Trazas/6 capas celofan en metal/ref metal celofan.txt')
# t_sam, E_sam = read_1file('./data/old/Paintmeter/Trazas/6 capas celofan en metal/sam metal celofan.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))
E_sam_max_idx = where(abs(E_sam) == E_sam_max)[0][0]


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))


ref_pulse_idx = centre_loc(E_ref)
k_bounds = [(1, 1e8),  # n1
            (1e-8, 1e8),  # k1
            (1e-8, 1e-3)  # thick1
            ]


res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam_w / E_ref_w, f_ref)
                             )

print(res)
# k = res.x
# pulses = list()
# for i in range(int(round(k.size / 2))):
#     pulses.append((k[2 * i + 1], k[2 * i]))
# pulses = sorted(pulses)
# delta_t = abs(pulses[0][0] - pulses[1][0])
# thickness = c_0 * delta_t * 1e-12 / (2 * 1.51 * cos(pi / 6))  # m
# thickness *= 1e3  # mm
# print(thickness, 'mm')
#
# plot(t_sam, E_sam, lw=1)
# E_sim = objective_function(k, E_ref, E_sam, t_ref)
# plot(t_ref, E_sim, lw=1)
# show()
