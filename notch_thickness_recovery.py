from TDSA import *
from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
from matplotlib.pyplot import *
from scipy.stats import norm
from scipy.signal import windows
from scipy.optimize import differential_evolution


# def thickness_val(n, order_m, f_ref, f_sam):
#     return c_0 * (order_m + 0.5) / (n - n_air) * (1 / f_ref - 1 / f_sam)
#
#
# def thickness_sigma2(n, sig_n, order_m, f_ref, sig_f_ref, f_sam, sig_f_sam):  # sigs in relative error (0 to 1)
#     aux1 = ((1 / f_ref - 1 / f_sam) * c_0 * (order_m + 0.5) / (n - n_air) ** 2) ** 2 * (sig_n * n) ** 2
#     aux2 = (c_0 * (order_m + 0.5) / (n - n_air)) ** 2 * ((sig_f_ref / f_ref) ** 2 + (sig_f_sam / f_sam) ** 2)
#     return aux1 + aux2
#
#
# def thickness_sigma2_aprox(n, order_m, f_ref, sig_f_ref, f_sam, sig_f_sam):  # sigs in relative error (0 to 1)
#     aux = (c_0 * (order_m + 0.5) / (n - n_air)) ** 2 * ((sig_f_ref / f_ref) ** 2 + (sig_f_sam / f_sam) ** 2)
#     return aux
#
#
# def objetive_func(x, loc, scale, a, m):
#     return cauchy.pdf(x, loc, scale) + a * x + m
#
#
# def cauchy_fit(params, *args):
#     loc, scale, a, m = params
#     x, y = args
#     return sum((cauchy.pdf(x, loc, scale) + a * x + m - y)**2)


def pdf_curve(freqs, loc, scale):
    # return cauchy.pdf(freqs, loc, scale)
    # return t.pdf(freqs, 10, loc, scale
    return exp(- abs(freqs - loc) / scale) / (2 * scale)


def H_sim(params, *args):
    # m2, m, b, loc1, scale1, loc2, scale2, loc3, scale3, loc4, scale4, loc5, scale5 = params
    # m2, m, b, loc1, scale1, loc2, scale2, loc3, scale3, loc4, scale4 = params
    # m2, m, b, loc1, scale1, loc2, scale2, loc3, scale3 = params
    # m2, m, b, loc1, scale1, loc2, scale2 = params
    m2, m, b, loc1, scale1 = params
    freqs, H_meas = args
    # base_line = m2 * freqs**2 + m / freqs + b
    base_line = m2 * freqs**2 + m * freqs + b
    # base_line = m * exp(m * freqs) + b
    # base_line = norm.pdf(freqs, 0, m) + b
    base_line -= pdf_curve(freqs, loc1, scale1)
    # base_line -= pdf_curve(freqs, loc2, scale2)
    # base_line -= pdf_curve(freqs, loc3, scale3)
    # base_line -= pdf_curve(freqs, loc4, scale4)
    # base_line -= pdf_curve(freqs, loc5, scale5)
    return base_line


def cost_function(params, *args):
    freqs, H_meas = args
    return sum((H_sim(params, *args) - H_meas)**2)


# t_ref, E_ref = read_slow_data('./data/notch/ref_notch.txt')
# t_sam, E_sam = read_slow_data('./data/notch/sam_notch.txt')
# t_ref, E_ref = read_1file('./data/marca_autodestructiva/ref.txt')
# t_ref, E_ref = read_1file('./data/marca_autodestructiva/edge.txt')
# t_sam, E_sam = read_1file('./data/marca_autodestructiva/sam.txt')
# t_sam, E_sam = read_1file('./data/marca_autodestructiva/edge.txt')
# t_ref, E_ref = read_slow_data('./data/notch/barrido 200 a 300 paso 0.1.txt')
# t_sam, E_sam = read_slow_data('./data/notch/barrido 300a 200 paso -0.1.txt')
# t_ref, E_ref = read_slow_data('./data/notch/new/ref1.txt')
# t_vid, E_vid = read_slow_data('./data/notch/new/vidre1.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper1.txt')
# t_ref, E_ref = read_slow_data('./data/notch/new/ref2.txt')
# t_vid, E_vid = read_slow_data('./data/notch/new/vidre2.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper2.txt')
t_ref, E_ref = read_slow_data('./data/notch/new/ref3.txt')
t_vid, E_vid = read_slow_data('./data/notch/new/vidre3.txt')
# t_pap, E_pap = read_slow_data('./data/notch/new/paper3.txt')
t_pap, E_pap = read_slow_data('./data/notch/new/celo0.txt')


win_idx = int(round(E_ref.size / 1.5))
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
w_filt = wiener_filter(E_ref_w)  # , beta=noise_floor(f_ref, E_ref_w, 1))
filt_size = int(E_ref_w.size / 6)
zero_size = E_ref_w.size - filt_size
tuk_filt = windows.tukey(filt_size, alpha=0.3)
tuk_filt = zero_padding(tuk_filt, 0, zero_size)
tuk_filt = roll(tuk_filt, 5)
# H_vid_w *= tuk_filt
# H_pap_w *= tuk_filt

f_min, f_max = f_min_max_idx(f_ref*1e12, 0.2, 0.4)
f_ref = f_ref[f_min:f_max]
H_vid_w = H_vid_w[f_min:f_max]
H_pap_w = H_pap_w[f_min:f_max]

# plot(f_ref, abs(H_vid_w), lw=1, label='vid')
# plot(f_ref, abs(H_pap_w), lw=1, label='pap')
# legend()
# show()
# quit()
# figure(1)
# # plot(t_ref, E_ref, lw=1)
# plot(t_vid, E_vid, lw=1)
# plot(t_pap, E_pap, lw=1)
#
bnds = list()
bnds.append((-10, 10))  # m2
bnds.append((-10, 10))  # m
bnds.append((-10, 10))  # b
bnds.append((f_ref[0], f_ref[-1]))  # loc1
# bnds.append((0.25, 0.35))  # loc1
bnds.append((0, 1))  # scale1
# bnds.append((f_ref[0], f_ref[-1]))  # loc2
# bnds.append((0.4, 0.55))  # loc2
# bnds.append((0, 1))  # scale2
# bnds.append((f_ref[0], f_ref[-1]))  # loc3
# bnds.append((0.6, 0.7))  # loc3
# bnds.append((0, 1))  # scale3
# bnds.append((f_ref[0], f_ref[-1]))  # loc4
# bnds.append((0, 0))  # scale4
# bnds.append((f_ref[0], f_ref[-1]))  # loc5
# bnds.append((0, 0))  # scale5
res_vid = differential_evolution(cost_function, bnds, args=(f_ref, toDb(H_vid_w)),
                             disp=True,
                             polish=True)
res_pap = differential_evolution(cost_function, bnds, args=(f_ref, toDb(H_pap_w)),
                             disp=True,
                             polish=True)
print(res_vid)
print(res_pap)

# m2_aux, m_aux, b_aux, loc_vid_1, scale1_aux, loc_vid_2, scale2_aux, loc_vid_3, scale3_aux, loc_vid_4, scale4_aux, loc_vid_5, scale5_aux= res_vid.x
# m2_aux, m_aux, b_aux, loc_pap_1, scale1_aux, loc_pap_2, scale2_aux, loc_pap_3, scale3_aux, loc_pap_4, scale4_aux, loc_pap_5, scale5_aux = res_pap.x
# locs_vid = sort(array([loc_vid_1, loc_vid_2, loc_vid_3, loc_vid_4, loc_vid_5]))
# locs_pap = sort(array([loc_pap_1, loc_pap_2, loc_pap_3, loc_pap_4, loc_pap_5]))
# m2_aux, m_aux, b_aux, loc_vid_1, scale1_aux, loc_vid_2, scale2_aux, loc_vid_3, scale3_aux, loc_vid_4, scale4_aux= res_vid.x
# m2_aux, m_aux, b_aux, loc_pap_1, scale1_aux, loc_pap_2, scale2_aux, loc_pap_3, scale3_aux, loc_pap_4, scale4_aux = res_pap.x
# locs_vid = sort(array([loc_vid_1, loc_vid_2, loc_vid_3, loc_vid_4]))
# locs_pap = sort(array([loc_pap_1, loc_pap_2, loc_pap_3, loc_pap_4]))
# m2_aux, m_aux, b_aux, loc_vid_1, scale1_aux, loc_vid_2, scale2_aux, loc_vid_3, scale3_aux = res_vid.x
# m2_aux, m_aux, b_aux, loc_pap_1, scale1_aux, loc_pap_2, scale2_aux, loc_pap_3, scale3_aux = res_pap.x
# locs_vid = sort(array([loc_vid_1, loc_vid_2, loc_vid_3])) * 1e12
# locs_pap = sort(array([loc_pap_1, loc_pap_2, loc_pap_3])) * 1e12
# m2_aux, m_aux, b_aux, loc_vid_1, scale1_aux, loc_vid_2, scale2_aux = res_vid.x
# m2_aux, m_aux, b_aux, loc_pap_1, scale1_aux, loc_pap_2, scale2_aux = res_pap.x
# locs_vid = sort(array([loc_vid_1, loc_vid_2]))
# locs_pap = sort(array([loc_pap_1, loc_pap_2]))
m2_aux, m_aux, b_aux, loc_vid_1, scale1_aux = res_vid.x
m2_aux, m_aux, b_aux, loc_pap_1, scale1_aux = res_pap.x
locs_vid = sort(array([loc_vid_1, loc_vid_1]))
locs_pap = sort(array([loc_pap_1, loc_pap_1]))


n_vid = 2.6
n_pap = 1.5
d_vid = 0.89*1e-3  # meters
m_orders = (d_vid * (n_vid + n_air) * locs_vid) / c_0
thicks = zeros(locs_vid.size)
for i in range(locs_vid.size):
    thicks[i] = ((c_0 * m_orders[i]) / (n_pap + n_air)) * (1 / locs_vid[i] - 1 / locs_pap[i])
    # thicks[i] = ((c_0 * (i + 1.5)) / (n_pap - n_air)) * (1 / locs_vid[i] - 1 / locs_pap[i])
    # thicks[i] = (- d_vid * (n_vid * locs_vid[i] - n_air * locs_pap[i]) * 1e12 + c_0 * m_orders[i]) / ((n_pap - n_air) * locs_pap[i] * 1e12)
# thicks *= 1e-12
print(thicks*1e6)
print(round(mean(thicks)*1e6), '+-', round(std(thicks)*1e6), 'um', '(' + str(round(100 * std(thicks) / mean(thicks))) + '%)')


figure(2)
plot(f_ref, toDb(H_vid_w), lw=1, label='vid')
plot(f_ref, H_sim(res_vid.x, f_ref, H_vid_w), lw=1, label='vft')
plot(f_ref, toDb(H_pap_w), lw=1, label='pap')
plot(f_ref, H_sim(res_pap.x, f_ref, H_pap_w), lw=1, label='pft')
xlim([f_ref[0], f_ref[-1]])
legend()


show()
