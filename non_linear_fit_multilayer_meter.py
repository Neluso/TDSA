from TDSA import *
from scipy.optimize import differential_evolution, NonlinearConstraint, LinearConstraint
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime
import math


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


reg_rel = 0.1
reg_matrix = identity(1 + layers * 4)
for i in range(reg_matrix.shape[0]):
    for j in range(reg_matrix.shape[1]):
        if i == j:
            reg_matrix[i, j] = 1
        elif i == 0:
            reg_matrix[i, j] = 0
        elif j == 0:
            reg_matrix[i, j] = 0
        else:
            reg_matrix[i, j] = reg_rel
reg_param = 0.1
reg_matrix *= reg_param



# debug variables
debug_value_1 = list()
debug_value_2 = list()
debug_value_3 = list()


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def ct2(n_l, n_l_1):
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1))
    n_l *= cos(theta(n_l))
    return - (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi) * exp(- k * phi)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def epsilon_cc(e_s, e_inf, tau, alpha, freq):  # Cole-Cole mode
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + (1j * omg * tau)**(alpha))
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def nk_from_eps_cc(e_s, e_inf, tau, alpha, freq):
    e_w = epsilon_cc(e_s, e_inf, tau, alpha, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(ns, ks, thicks, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, ns[-1] - 1j * ks[-1]) * ones(freq.size)
    # # real indices
    # rlm1l = cr_l_1_l(ns[-1], ns[-2])
    # tt = ct2(ns[-1] , ns[-2])
    # complex indices
    rlm1l = cr_l_1_l(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    tt = ct2(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    exp_phi = phase_factor(ns[-1], ks[-1], thicks[-1], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # inner layers, no contact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    for i in range(1, layers - 1):
        i = layers - i - 1
        # # real indices
        # rlm1l = cr_l_1_l(ns[i], ns[i - 1])
        # tt = ct2(ns[i], ns[i - 1])
        # complex indices
        rlm1l = cr_l_1_l(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        tt = ct2(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
        H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # most outer layer, in contact with air
    # real indices
    # rlm1l = cr_l_1_l(ns[0], n_air)
    # tt = ct2(ns[0], n_air)
    # complex indices
    rlm1l = cr_l_1_l(ns[0] - 1j * ks[0], n_air_cplx)
    tt = ct2(ns[0] - 1j * ks[0], n_air_cplx)
    exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


def cost_function(params, *args):
    ns = list()
    ks = list()
    thicks = list()
    H_meas, freqs, E_ref_w = args
    for i in range(layers):
        e_s = params[4 * i]
        e_inf = params[4 * i + 1]
        tau = params[4 * i + 2]
        thicks.append(params[4 * i + 3])
        n, k = nk_from_eps(e_s, e_inf, tau, freqs)  # debye model
        # n, k = nk_from_eps_cc(e_s, e_inf, tau, 0.999, freqs)  # Cole-Cole model
        ns.append(n)
        ks.append(k)
    H_teo = H_sim(array(ns), array(ks), array(thicks), freqs)
    logs = log(abs(H_teo)) - log(abs(H_meas))
    angles = abs(unwrap(angle(H_teo))) - abs(unwrap(angle(H_meas)))
    logs = abs(logs)
    angles = abs(angles)

    reg_params = dot(reg_matrix, params)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)

    reg_vals = E_sam_teo - E_sam

    return dot(reg_vals, reg_vals) + dot(reg_params, reg_params)
    # return sum(abs(logs + 0.1 * angles))


# def cons(params, *args):
#     thicks = list()
#     for i in range(layers):
#         thicks.append(params[3 * i + 2])
#     thicks = array(thicks)
#     # H_meas, freqs = args
#     return array([thicks[0] - thicks[1], thicks[1] - thicks[2]])


# Main script
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')

print('DSP')
delta_t_ref = mean(diff(t_ref))
ref_pulse_idx = centre_loc(E_ref)
window = tukey(2 * ref_pulse_idx)
window = zero_padding(window, 0, E_ref.size - window.size)
E_ref *= window
E_sam *= window
enlargement = 2 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))
f_min, f_max = f_min_max_idx(f_ref, 0.2, 0.6)
# f_min, f_max = f_min_max_idx(f_ref, 0, f_ref[-1]*1e-12)
# f_min, f_max = f_min_max_idx(f_ref, 0, 10)
f_ref = f_ref[f_min:f_max]
E_ref_w = E_ref_w[f_min:f_max]
f_sam = f_sam[f_min:f_max]
E_sam_w = E_sam_w[f_min:f_max]

H_w = E_sam_w / E_ref_w
filt = wiener_filter(E_ref_w)
H_w *= filt


cons = array(  # constraint matrix
    [
        [
            0, 0, 1, 0, 0, -1, 0, 0, 0
        ],
        [
            0, 0, 0, 0, 0, 1, 0, 0, -1
        ]
    ]
)


fh = open('./output/epsilons/log.txt', 'a')

print('Fitting')
alpha = 1  # 7
beta = 50  # 15
k_bounds = [  # blind calibration
    (0, 100e-6),  # air thickness
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6)
]
# k_bounds = [  # blind calibration
#     (alpha, beta), (alpha, beta), (0, 1), (1e-14, 1e-12), (0, 100e-6),
#     (alpha, beta), (alpha, beta), (0, 1), (1e-14, 1e-12), (0, 100e-6),
#     (alpha, beta), (alpha, beta), (0, 1), (1e-14, 1e-12), (0, 100e-6)
# ]
# k_bounds = [  # calibration
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (81e-6, 81e-6),  # 87*1e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (19e-6, 19e-6),  # 20*1e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (5e-6, 5e-6)  # 5*1e-6)
# ]

t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(H_w, f_ref, E_ref_w),
                             # strategy='best1bin',
                             popsize=90,
                             # maxiter=2000,
                             disp=True,  # step-by-step cost_function value
                             # mutation=1.5,
                             # workers=2,
                             # constraints=NonlinearConstraint(cons,0,ones(2)),
                             # constraints=(LinearConstraint(cons, 1e-12*ones(2), ones(2))),
                             polish=True
                             )


t2 = time_ns()
print('alpha =', alpha)
print('beta =', beta)
print(res)
print()
print('Results:')
n_w, k_w = nk_from_eps(res.x[0], res.x[1], res.x[2], f_ref)
n_g, k_g = nk_from_eps(res.x[4], res.x[5], res.x[6], f_ref)
d_w = res.x[3]
d_g = res.x[7]
print('White coat -', 'n:', round(mean(n_w), 2), 'k:', round(mean(k_w), 2), 'd:', round(d_w * 1e6, 2), 'um')
print('Green coat -', 'n:', round(mean(n_g), 2), 'k:', round(mean(k_g), 2), 'd:', round(d_g * 1e6, 2), 'um')
if layers == 3:
    n_p, k_p = nk_from_eps(res.x[8], res.x[9], res.x[10], f_ref)
    d_p = res.x[11]
    print('Primer coat -', 'n:', round(mean(n_p), 2), 'k:', round(mean(k_p), 2), 'd:', round(d_p * 1e6, 2), 'um')
    print('Total:', round((d_w + d_g + d_p) * 1e6, 2), 'um')
    ns = array([n_w, n_g, n_p])
    ks = array([k_w, k_g, k_p])
    thikcs = array([d_w, d_g, d_p])
else:
    print('Total:', round((d_w + d_g) * 1e6, 2), 'um')
    ns = array([n_w, n_g])
    ks = array([k_w, k_g])
    thikcs = array([d_w, d_g])
print()
secs = (t2 - t1) * 1e-9
mins_aux = secs / 60
mins = int(trunc(mins_aux))
secs = ((mins_aux - mins) * 60)
secs = int(round(secs, 0))
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()
fh.write('alpha = ' + str(alpha) + '\n')
fh.write('beta = ' + str(beta) + '\n')
fh.write(str(res) + '\n')
fh.write('Results:\n')
n_w, k_w = nk_from_eps(res.x[0], res.x[1], res.x[2], f_ref)
n_g, k_g = nk_from_eps(res.x[4], res.x[5], res.x[6], f_ref)
d_w = res.x[3]
d_g = res.x[7]
fh.write('White coat -\t' + 'n: ' + str(round(mean(n_w), 2)) + '\tk: ' + str(round(mean(k_w), 2)) + '\td: ' + str(round(d_w * 1e6, 2)) + ' um' + '\n')
fh.write('Green coat -\t' + 'n: ' + str(round(mean(n_g), 2)) + '\tk: ' + str(round(mean(k_g), 2)) + '\td: ' + str(round(d_g * 1e6, 2)) + ' um' + '\n')
n_p, k_p = nk_from_eps(res.x[8], res.x[9], res.x[10], f_ref)
d_p = res.x[11]
fh.write('Primer coat -\t' + 'n: ' + str(round(mean(n_p), 2)) + '\tk:' + str(round(mean(k_p), 2)) + '\td:' + str(round(d_p * 1e6, 2)) + ' um' + '\n')
fh.write('Total: ' + str(round((d_w + d_g + d_p) * 1e6, 2)) + ' um' + '\n')
ns = array([n_w, n_g, n_p])
ks = array([k_w, k_g, k_p])
thikcs = array([d_w, d_g, d_p])
fh.write('\n')
secs = (t2 - t1) * 1e-9
mins_aux = secs / 60
mins = int(trunc(mins_aux))
secs = ((mins_aux - mins) * 60)
secs = int(round(secs, 0))
fh.write('Processing time (mm:ss): ' + str(mins) + ':' + str(secs) + '\n\n')


# Plotting
H_teo = H_sim(ns, ks, thikcs, f_ref)
f_ref *= 1e-12
f_sam *= 1e-12

fig, axs = subplots(2)
fig.suptitle('Abs/Phase')

axs[0].plot(f_ref, abs(H_w), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_w)), lw=1)

axs[0].plot(f_ref, abs(H_teo), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_teo)), lw=1)


axs[0].set_ylabel(r'$\rho$')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], f_ref[-1]])
axs[1].set_ylabel(r'$\phi \ (rad)$')
axs[1].set_xlim([f_ref[0], f_ref[-1]])
xlabel(r'$f\ (THz)$')
savefig('./output/epsilons/abs_phi_a'+str(alpha)+'_b'+str(beta)+'.png', format='PNG')

fig, axs = subplots(2)
fig.suptitle('Re/Im')

axs[0].plot(f_ref, real(H_w), lw=1)
axs[1].plot(f_ref, imag(H_w), lw=1)

axs[0].plot(f_ref, real(H_teo), lw=1)
axs[1].plot(f_ref, imag(H_teo), lw=1)


axs[0].set_ylabel(r'$Re$')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], f_ref[-1]])
axs[1].set_ylabel(r'$Im$')
axs[1].set_xlim([f_ref[0], f_ref[-1]])
xlabel(r'$f\ (THz)$')
savefig('./output/epsilons/re_im_a' + str(alpha) + '_b' + str(beta) + '.png', format='PNG')

figure(3)
plot(f_ref, n_w, lw=1, label='White')
plot(f_ref, n_g, lw=1, label='Green')
if layers == 3:
    plot(f_ref, n_p, lw=1, label='Primer')
legend()
savefig('./output/epsilons/ns_a' + str(alpha) + '_b' + str(beta) + '.png', format='PNG')
# clf()

show()
