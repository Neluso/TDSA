# We implement here a simple modification for the H_sim function, proposed by Mechelen et al.
# http://dx.doi.org/10.1364/OL.39.003853


from TDSA import *
from scipy.optimize import differential_evolution, minimize
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime
import datetime
import math


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


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


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k

def H_sim(ns, ks, thicks, d_air, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, ns[-1] - 1j * ks[-1]) * ones(freq.size)
    rlm1l = cr_l_1_l(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    tt = ct2(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    exp_phi = phase_factor(ns[-1], ks[-1], thicks[-1], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # inner layers, no contact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    for i in range(1, layers - 1):
        i = layers - i - 1
        rlm1l = cr_l_1_l(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        tt = ct2(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
        H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # most outer layer, in contact with air
    rlm1l = cr_l_1_l(ns[0] - 1j * ks[0], n_air_cplx)
    tt = ct2(ns[0] - 1j * ks[0], n_air_cplx)
    exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return phase_factor(n_air, 0, d_air, freq) * H_i


def cost_function(params, *args):
    d_air = params[0]
    params = params[1:]
    ns = list()
    ks = list()
    thicks = list()
    E_sam, E_ref_w, freqs = args
    for i in range(layers):
        e_s = params[4 * i]
        e_inf = params[4 * i + 1]
        tau = params[4 * i + 2]
        thicks.append(params[4 * i + 3])
        n, k = nk_from_eps(e_s, e_inf, tau, freqs)  # debye model
        ns.append(n)
        ks.append(k)
    ns = array(ns)
    ks = array(ks)
    thicks = array(thicks)
    H_teo = H_sim(ns, ks, thicks, d_air, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    return sum((E_sam_teo - E_sam)**2)


# Main script
# Boleto 176054
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat1_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/sam metal wcoat 2_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')


print('DSP')
delta_t_ref = mean(diff(t_ref))
ref_pulse_idx = centre_loc(E_ref)
window = tukey(E_ref.size)
E_ref *= window
E_sam *= window
t_ref *= 1e-12
t_sam *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
H_w = E_sam_w / E_ref_w

alpha = 7
beta = 12
# k_bounds = [  # calibration
#     (0, 100e-6),  # air thickness
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6)
# ]
k_bounds = [  # calibration
    (0, 100e-6),  # air thickness
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (81e-6, 81e-6),  # 87*1e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (19e-6, 19e-6),  # 20*1e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (5e-6, 5e-6)  # 5*1e-6)
]

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

# TODO review full H_sim
# i = 0
# strategies = ['best1bin', 'best1exp', 'rand1exp', 'randtobest1exp', 'currenttobest1exp', 'best2exp', 'rand2exp',
#               'randtobest1bin', 'currenttobest1bin', 'best2bin', 'rand2bin', 'rand1bin']
print('Fitting')
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             strategy='best1bin',
                             # popsize=150,
                             # maxiter=2000,
                             disp=True,  # step cost_function value
                             # mutation=1.5,
                             workers=1,
                             # constraints=NonlinearConstraint(cons,,ones(2)),
                             # constraints=LinearConstraint(cons, 1e-12*ones(2), ones(2)),
                             polish=False
                             )
t2 = time_ns()
print()
d_air = res.x[0]
print(res)
res.x = res.x[1:]
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

H_teo = H_sim(ns, ks, thikcs, d_air, f_ref)
E_sam_teo = irfft(H_teo * E_ref_w, n=t_sam.size)


f_ref *= 1e-12
f_sam *= 1e-12
t_ref *= 1e12
t_sam *= 1e12

fig, axs = subplots(2)
fig.suptitle('Abs/Phase')

axs[0].plot(f_ref, abs(H_w), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_w)), lw=1)

axs[0].plot(f_ref, abs(H_teo), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_teo)), lw=1)


axs[0].set_ylabel(r'$\rho$')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], 1.2])
axs[0].set_ylim([0, 1.3])
axs[1].set_ylabel(r'$\phi \ (rad)$')
axs[1].set_xlim([f_ref[0], 1.2])
axs[1].set_ylim([- 4 * pi, pi])
xlabel(r'$f\ (THz)$')


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


figure(5)
plot(t_sam, E_sam, lw=1)
plot(t_sam, E_sam_teo, lw=1)
figure(6)
plot(t_sam, E_sam - E_sam_teo, lw=1)


show()
