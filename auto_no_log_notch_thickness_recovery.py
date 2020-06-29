from TDSA import *
from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
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
    return sqrt(0.25 * ((ct_total**2 + 1) + 2*(ct_total**2 - 1)*erf_func + (ct_total**2 + 1)*erf_func**2 + 2*(1 - erf_func**2)*cos_func))


def ct(n_l, n_l_1):
    return 2 * n_l / (n_l + n_l_1)


def abs_E_ref_w(erf_par, d_ref, freqs, n_ref, k_ref, bias):
    phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air)) / c_0
    erf_func = erf(erf_par * freqs)
    cos_func = exp(- 2 * pi * freqs * k_ref * d_ref / c_0) * cos(phi12)
    ct_tot = ct(n_ref, n_air) * ct(n_air, n_ref)
    # cos_func = cos(phi12)
    # this is not even my
    return final_form(erf_func, cos_func, ct_tot) - bias


def abs_E_obj_w(erf_par, d_obj, freqs, d_ref, n_ref, k_ref, n_obj, bias):
    phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air) + d_obj * (n_obj - n_air)) / c_0  # sobre el vidre
    # phi12 = 2 * pi * freqs * (d_ref * (n_ref - n_air) - d_obj * (n_obj - n_air)) / c_0  # contra el vidre
    erf_func = erf(erf_par * freqs)
    cos_func = exp(- 2 * pi * freqs * k_ref * d_ref / c_0) * cos(phi12)
    ct_tot = ct(n_air, n_ref) * ct(n_ref, n_obj) * ct(n_obj, n_air)
    # ct_tot = ct(n_air, n_ref) * ct(n_obj, n_air)
    # cos_func = cos(phi12)
    return final_form(erf_func, cos_func, ct_tot) - bias


def cost_function_ref(params, *args):
    erf_par, d_ref, bias = params
    freqs, H_meas_ref, n_ref, k_ref = args
    abs_ref = abs_E_ref_w(erf_par, d_ref, freqs, n_ref, k_ref, bias)
    return sum((abs_ref - H_meas_ref)**2) + 0.001 * sum(abs(params)**2)


def cost_function_obj(params, *args):
    erf_par, d_obj, n_obj, bias = params
    freqs, H_meas_obj, d_ref, n_ref, k_ref = args
    abs_obj = abs_E_obj_w(erf_par, d_obj, freqs, d_ref, n_ref, k_ref, n_obj, bias)
    return sum((abs_obj - H_meas_obj)**2) + 0.001 * sum(abs(params)**2)


bnds_ref = list()
bnds_ref.append((0, 1e-11))  # erf_par_vid
bnds_ref.append((0.7e-3, 1.2e-3))  # d_ref
bnds_ref.append((0, 1))  # bias

bnds_obj = list()
bnds_obj.append((0, 1e-11))  # erf_par_pap
bnds_obj.append((1e-6, 1000e-6))  # d_obj
bnds_obj.append((1.35, 1.8))  # n_obj
bnds_obj.append((0, 1))  # bias

# n_paps = [1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39]
# n_paps = [1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.51]
# n_paps = [1.52, 1.53, 1.54, 1.55, 1.56]

wh = open('./output/notch_results.csv', 'w')

for i in range(10):
    t_ref, E_ref = read_slow_data('./data/notch/camp1/ref'+str(i+1)+'.txt')
    t_vid, E_vid = read_slow_data('./data/notch/camp1/vidre'+str(i+1)+'.txt')
    t_pap, E_pap = read_slow_data('./data/notch/camp1/celo'+str(i+1)+'.txt')
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
    
    f_min, f_max = f_min_max_idx(f_ref * 1e12, 0.2, 1.0)
    f_ref = f_ref[f_min:f_max]
    H_vid_w = H_vid_w[f_min:f_max]
    H_pap_w = H_pap_w[f_min:f_max]

    n_vid_itpl, k_vid_itpl = read_calibrations()
    n_vid = n_vid_itpl(f_ref * 1e12)
    k_vid = k_vid_itpl(f_ref * 1e12)
    
    res_ref = differential_evolution(cost_function_ref, bnds_ref, args=(f_ref * 1e12, abs(H_vid_w), n_vid, k_vid),
                                     # disp=True,
                                     # popsize=600,
                                     # maxiter=5000,
                                     polish=True)
    d_vid = res_ref.x[1]
    res_obj = differential_evolution(cost_function_obj, bnds_obj,
                                     args=(f_ref * 1e12, abs(H_pap_w), d_vid, n_vid, k_vid),
                                     # disp=True,
                                     # popsize=600,
                                     # maxiter=5000,
                                     polish=True)
    d_pap = res_obj.x[1]
    n_pap = res_obj.x[2]
    wh.write(str(round(n_pap, 2)) + ';' + str(round(d_pap * 1e6, 0)).split('.')[0] + '\n')
    print('Sample', i+1)
    print('d_vid =', round(d_vid * 1e3, 2), 'mm', '=', round(d_vid * 1e6, 0), 'um')
    print('d_obj =', round(d_pap * 1e6, 0), 'um')
    print('n_obj =', n_pap)
    print(res_obj.x)
    print()
