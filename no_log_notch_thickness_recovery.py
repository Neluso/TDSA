from TDSA import *
from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
from matplotlib.pyplot import *
from scipy.signal import windows
from scipy.optimize import differential_evolution
from scipy.special import erf
from scipy.interpolate import Akima1DInterpolator


def read_calibrations():
    file_path = './data/calibrations/nvidrio.txt'
    fh = open(file_path, 'r')
    lines = fh.read().split('\n')
    freqs1 = zeros(len(lines) - 1)
    vals1 = zeros(len(lines) - 1)
    for i in range(len(lines)):
        line = lines[i]
        line = line.split(',')
        try:
            freqs1[i] = float(line[0])
            vals1[i] = float(line[1])
        except:
            continue
    fh.close()
    file_path = './data/calibrations/kvidrio.txt'
    fh = open(file_path, 'r')
    lines = fh.read().split('\n')
    freqs2 = zeros(len(lines) - 1)
    vals2 = zeros(len(lines) - 1)
    for i in range(len(lines)):
        line = lines[i]
        line = line.split(',')
        try:
            freqs2[i] = float(line[0])
            vals2[i] = float(line[1])
        except:
            continue
    return Akima1DInterpolator(freqs1, vals1), Akima1DInterpolator(freqs2, vals2),


# def final_form(erf_func, cos_func):
#     # hem obtés: sqrt(0.5 * (1 + erf_func + (1 - erf_func) * cos_func))
#     # return sqrt(0.5 * (1 - erf_func + (1 - erf_func) * cos_func))  # funciona
#     return sqrt(0.5 * (1 + erf_func + (1 - erf_func) * cos_func))


def final_form(erf_func, cos_func, ct_total):
    # hem obtés: sqrt(0.5 * (1 + erf_func + (1 - erf_func) * cos_func))
    # return sqrt(0.5 * (1 - erf_func + (1 - erf_func) * cos_func))  # funciona
    return sqrt(0.25 * ((abs(ct_total)**2 + 1) + 2*(abs(ct_total)**2 - 1)*erf_func + (abs(ct_total)**2 + 1)*erf_func**2 + 2*(1 - erf_func**2)*cos_func))


def ct(n_l, n_l_1):
    return 2 * n_l / (n_l + n_l_1)


def abs_E_ref_w(erf_par, d_ref, freqs, n_ref, k_ref, bias):
    phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air)) / c_0
    erf_func = erf(erf_par * freqs)
    ct_tot = ct(n_ref - 1j * k_ref * freqs, n_air)**2  # * ct(n_air, n_ref - 1j * k_ref * freqs)
    phi12 += angle(ct_tot)
    cos_func = exp(- 2 * pi * freqs * k_ref * d_ref / c_0) * cos(phi12)
    # cos_func = cos(phi12)
    # this is not even my
    return final_form(erf_func, cos_func, ct_tot) - bias


def abs_E_obj_w(erf_par, d_obj, freqs, d_ref, n_ref, k_ref, n_obj, k_obj, bias):
    phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air) + d_obj * (n_obj - n_air)) / c_0  # sobre el vidre
    # phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air) - d_obj * (n_obj - n_air)) / c_0  # contra el vidre
    erf_func = erf(erf_par * freqs)
    ct_tot = ct(n_air, n_ref - 1j * k_ref)
    ct_tot = ct_tot * ct(n_ref - 1j * k_ref * freqs, n_obj - 1j * k_obj * freqs)
    ct_tot = ct_tot * ct(n_obj - 1j * k_obj * freqs, n_air)
    phi12 += angle(ct_tot)
    cos_func = exp(- 2 * pi * freqs * k_ref * d_ref / c_0) * cos(phi12)
    # cos_func = cos(phi12)
    return final_form(erf_func, cos_func, ct_tot) - bias


def cost_function_ref(params, *args):
    erf_par, d_ref, bias = params
    freqs, H_meas_ref, n_ref, k_ref = args
    abs_ref = abs_E_ref_w(erf_par, d_ref, freqs, n_ref, k_ref, bias)
    return sum((abs_ref - H_meas_ref)**2) + 0.01 * sum(abs(params)**2)


def cost_function_obj(params, *args):
    erf_par, d_obj, n_obj, k_obj, bias = params
    freqs, H_meas_obj, d_ref, n_ref, k_ref = args
    abs_obj = abs_E_obj_w(erf_par, d_obj, freqs, d_ref, n_ref, k_ref, n_obj, k_obj, bias)
    return sum((abs_obj - H_meas_obj)**2) + 0.01 * sum(abs(params)**2)


# t_ref, E_ref = read_slow_data('./data/notch/new/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/new/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper1.txt')
# t_ref, E_ref = read_slow_data('./data/notch/new/ref2.txt')
# t_vid, E_vid = read_slow_data('./data/notch/new/vidre2.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper2.txt')
# t_ref, E_ref = read_slow_data('./data/notch/new/ref3.txt')
# t_vid, E_vid = read_slow_data('./data/notch/new/vidre3.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper3.txt')

# t_pap, E_pap = read_slow_data('./data/notch/new/celo0.txt')


# t_ref_calib, E_ref_calib = read_1file('./data/notch/notch_celo/ref_spec.txt')
# t_sam_calib, E_sam_calib = read_1file('./data/notch/notch_celo/vidre_spec.txt')
# 1 celo
# t_ref, E_ref = read_slow_data('./data/notch/notch_celo/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/notch_celo/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/notch_celo/celo1.txt')
# t_ref, E_ref = read_slow_data('./data/notch/notch_celo/ref2.txt')
# t_vid, E_vid = read_slow_data('./data/notch/notch_celo/vidre2.txt')
# t_pap, E_pap = read_slow_data('./data/notch/notch_celo/celo2.txt')
# t_ref, E_ref = read_slow_data('./data/notch/notch_celo/ref3.txt')
# t_vid, E_vid = read_slow_data('./data/notch/notch_celo/vidre3.txt')
# t_pap, E_pap = read_slow_data('./data/notch/notch_celo/celo3.txt')
# 2 celos
# t_ref, E_ref = read_slow_data('./data/notch/2celo/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/2celo/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/2celo/celo1.txt')
# 3 celos
# t_ref, E_ref = read_slow_data('./data/notch/3celo/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/3celo/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/3celo/celo1.txt')
# 4 celos
# t_ref, E_ref = read_slow_data('./data/notch/4celo/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/4celo/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/4celo/celo1.txt')

# Campanya 1 (2 celos)
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo1.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref2.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre2.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo2.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref3.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre3.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo3.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref4.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre4.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo4.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref5.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre5.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo5.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref6.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre6.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo6.txt')
t_ref, E_ref = read_slow_data('./data/notch/camp1/ref7.txt')
t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre7.txt')
t_pap, E_pap = read_slow_data('./data/notch/camp1/celo7.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref8.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre8.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo8.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref9.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre9.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo9.txt')
# t_ref, E_ref = read_slow_data('./data/notch/camp1/ref10.txt')
# t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre10.txt')
# t_pap, E_pap = read_slow_data('./data/notch/camp1/celo10.txt')

plot(t_vid, E_vid)
show()
quit()


win_idx = int(round(E_ref.size / 3))
window = windows.tukey(win_idx, 0.2)
window = zero_padding(window, 0, E_ref.size - win_idx)


E_ref *= window
E_vid *= window
E_pap *= window


delta_t_ref = mean(diff(t_ref))
enlargement = 2 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_vid = zero_padding(E_vid, 0, enlargement)
t_vid = concatenate((t_vid, t_vid[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_pap = zero_padding(E_pap, 0, enlargement)
t_pap = concatenate((t_pap, t_pap[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_vid, E_vid_w = fourier_analysis(t_vid, E_vid)
f_pap, E_pap_w = fourier_analysis(t_pap, E_pap)
H_vid_w = E_vid_w / E_ref_w
H_pap_w = E_pap_w / E_ref_w
w_filt = wiener_filter(E_ref_w, beta=1e-8)  # , beta=noise_floor(f_ref, E_ref_w, 1))
filt_size = int(E_ref_w.size / 6)
zero_size = E_ref_w.size - filt_size
tuk_filt = windows.tukey(filt_size, alpha=0.3)
tuk_filt = zero_padding(tuk_filt, 0, zero_size)
tuk_filt = roll(tuk_filt, 30)
# H_vid_w *= w_filt
# H_pap_w *= w_filt


f_min, f_max = f_min_max_idx(f_ref*1e12, 0.15, 1.1)
f_ref = f_ref[f_min:f_max]
H_vid_w = H_vid_w[f_min:f_max]
H_pap_w = H_pap_w[f_min:f_max]
E_pap_w = E_pap_w[f_min:f_max]


bnds_ref = list()
bnds_ref.append((0, 1e-11))  # erf_par_vid
bnds_ref.append((0.7e-3, 1.2e-3))  # d_ref
# bnds_ref.append((0.95e-3, 0.95e-3))  # d_ref
bnds_ref.append((0, 1))  # bias

bnds_obj = list()
bnds_obj.append((0, 1e-11))  # erf_par_pap
bnds_obj.append((1e-6, 1000e-6))  # d_obj
bnds_obj.append((1.35, 1.8))  # n_obj
bnds_obj.append((0, 1))  # k_obj
bnds_obj.append((0, 1))  # bias

n_vid_itpl, k_vid_itpl = read_calibrations()
n_vid = n_vid_itpl(f_ref*1e12)
k_vid = k_vid_itpl(f_ref*1e12)
# n_vid = 2.6  # * (1 - 0.01 * f_ref)
# k_vid = 0.0  # * f_ref
# n_pap = 1.51 - 0.02 * f_ref
# n_pap = 1.42  # indice ajustado a mano...

res_ref = differential_evolution(cost_function_ref, bnds_ref, args=(f_ref*1e12, abs(H_vid_w), n_vid, k_vid),
                                 disp=True,
                                 # popsize=600,
                                 # maxiter=5000,
                                 polish=True)
d_vid = res_ref.x[1]

figure(2)
title('Vidre')
plot(f_ref, abs(H_vid_w), lw=1, label='meas')
plot(f_ref, abs_E_ref_w(res_ref.x[0], d_vid, f_ref*1e12, n_vid, k_vid, res_ref.x[2]), lw=1, label='fit')
xlim([f_ref[0], f_ref[-1]])
ylim([0, 1])
legend()

# print(res_ref)
# print('d_vid =', round(d_vid*1e3, 2), 'mm', '=', round(d_vid*1e6, 0), 'um')
# show()
# quit()

res_obj = differential_evolution(cost_function_obj, bnds_obj, args=(f_ref*1e12, abs(H_pap_w), d_vid, n_vid, k_vid),
                                 disp=True,
                                 # popsize=600,
                                 # maxiter=5000,
                                 polish=True)

d_pap = res_obj.x[1]
n_pap = res_obj.x[2]

print(res_ref)
print(res_obj)

print('d_vid =', round(d_vid*1e3, 2), 'mm', '=', round(d_vid*1e6, 0), 'um')
print('d_pap =', round(d_pap*1e6, 0), 'um')
print('n_pap =', n_pap)


figure(3)
title('Vidre+Celo')
plot(f_ref, abs(H_pap_w), lw=1, label='meas')
plot(f_ref, abs_E_obj_w(res_obj.x[0], d_pap, f_ref*1e12, d_vid, n_vid, k_vid, n_pap, res_obj.x[3], res_obj.x[4]), lw=1, label='fit')
xlim([f_ref[0], f_ref[-1]])
ylim([0, 1])
legend()


show()
