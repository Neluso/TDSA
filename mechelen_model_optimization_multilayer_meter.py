# We implement here a simple modification for the H_sim function, proposed by Mechelen et al.
# http://dx.doi.org/10.1364/OL.39.003853


from TDSA import *
from scipy.optimize import differential_evolution
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0
# n_subs = 3.00 - 1.0 * 1j  # substrate refractive index -- cork 3.0


reg_rel = 0.1
reg_matrix = identity(1 + layers * 3)
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
print(reg_matrix.shape)
# quit()


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


def cost_function(k, *args):
    E_sam, E_ref_w, freqs = args
    reg_params = dot(reg_matrix, k)
    d_air = k[0]
    k = k[1:]
    ns = list()
    ks = list()
    thicks = list()
    for i in range(layers):
        ns.append(k[3*i])
        ks.append(k[3*i + 1])
        thicks.append(k[3*i + 2])
    ns = array(ns)
    ks = array(ks)
    thicks = array(thicks)
    
    H_teo = H_sim(ns, ks, thicks, d_air, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    
    reg_vals = E_sam_teo - E_sam
    
    return dot(reg_vals, reg_vals) + dot(reg_params, reg_params)


# def cons(params, *args):
#     thicks = list()
#     for i in range(layers):
#         thicks.append(params[3 * i + 2])
#     thicks = array(thicks)
#     # H_meas, freqs = args
#     return array([thicks[0] - thicks[1], thicks[1] - thicks[2]])


# Main script
# Boleto 176054
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat2_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')


print('DSP')
delta_t_ref = mean(diff(t_ref))
ref_pulse_idx = centre_loc(E_ref)
window = tukey(2 * ref_pulse_idx)
window = zero_padding(window, 0, E_ref.size - window.size)
E_ref *= window
E_sam *= window
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
# E_sam_w *= wiener_filter(E_ref_w)

delta_f_ref = mean(diff(f_ref))
# f_min, f_max = f_min_max_idx(f_ref, 0.1, 1.1)
# f_ref = f_ref[f_min:f_max]
# E_ref_w = E_ref_w[f_min:f_max]
# f_sam = f_sam[f_min:f_max]
# E_sam_w = E_sam_w[f_min:f_max]

H_w = E_sam_w / E_ref_w
# H_w /= amax(H_w)


# k_bounds = [  # narrow
#     (0, 100e-6),  # air thickness
#     (1.8, 5), (0, 10), (50*1e-6, 150*1e-6),
#     (1.8, 5), (0, 10), (5*1e-6, 50*1e-6),
#     (1.8, 5), (0, 10), (1*1e-6, 10*1e-6)
# ]
# k_bounds = [  # tailored optical parameters (emulating a calibration)
#     (0, 100e-6),  # air thickness
#     (3, 4), (0, 1), (0, 1e-3),
#     (3.5, 4.1), (6.9, 8.1), (0, 1e-3),
#     (4, 4.4), (7, 10), (0, 1e-3)
# ]
# k_bounds = [  # highly tailored optical parameters
#     (1e-6, 50e-6),  # air thickness
#     (3, 4), (0, 1), (20e-6, 50e-6),
#     (3.5, 4.1), (6.9, 8.1), (40e-6, 53e-6),
#     (4, 4.4), (7, 10), (5e-6, 10e-6)
# ]
# k_bounds = [  # very narrow
#     (0, 100e-6),  # air thickness
#     (3, 4), (0.1, 0.4), (50*1e-6, 150*1e-6),
#     (3.5, 4.1), (6.9, 8.1), (5*1e-6, 50*1e-6),
#     (4, 4.4), (7, 10), (1*1e-6, 10*1e-6)
# ]
k_bounds = [  # wide
    (0, 100e-6),  # air thickness
    (1.2, 10), (0, 10), (0, 1e-3),
    (1.2, 10), (0, 10), (0, 1e-3),
    (1.2, 10), (0, 10), (0, 1e-3)
]
# k_bounds = [  # testing
#     (0, 100e-6),  # air thickness
#     (1, 100), (0, 10), (81e-6, 81e-6),
#     (1, 100), (0, 10), (19e-6, 19e-6),
#     (1, 100), (0, 10), (4e-6, 4e-6)
# ]
# k_bounds = [  # calibrating 177910
#     (0, 100e-6),  # air thickness
#     (1, 100), (0, 10), (34e-6, 34e-6),
#     (1, 100), (0, 10), (47e-6, 47e-6),
#     (1, 100), (0, 10), (7e-6, 7e-6)
# ]
# k_bounds = [  # reconstructicting 177910
#     (0, 100e-6),  # air thickness
#     (3.5, 4), (0, 1), (30e-6, 40e-6),
#     (1.5, 3), (0, 5), (40e-6, 50e-6),
#     (20, 30), (0, 5), (5e-6, 10e-6)
# ]
# k_bounds = [  # calibrating 180881
#     (0, 100e-6),  # air thickness
#     (1, 100), (0, 10), (47e-6, 47e-6),
#     (1, 100), (0, 10), (27e-6, 27e-6),
#     (1, 100), (0, 10), (11e-6, 11e-6)
# ]
# k_bounds = [  # reconstructicting 180881
#     (0, 100e-6),  # air thickness
#     (3, 4), (0, 1), (40e-6, 50e-6),
#     (2, 4), (0, 10), (30e-6, 40e-6),
#     (10, 20), (0, 5), (5e-6, 15e-6)
# ]

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
                             # popsize=150,
                             # maxiter=2000,
                             disp=True,  # step cost_function value
                             polish=True
                             )
t2 = time_ns()
print()
d_air = res.x[0]
print(res)
res.x = res.x[1:]
print('Results:')
print('White coat -', 'n:', round(res.x[0], 2), 'k:', round(res.x[1], 2), 'd:', round(res.x[2] * 1e6, 2), 'um')
print('Green coat -', 'n:', round(res.x[3], 2), 'k:', round(res.x[4], 2), 'd:', round(res.x[5] * 1e6, 2), 'um')
if layers == 3:
    print('Primer coat -', 'n:', round(res.x[6], 2), 'k:', round(res.x[7], 2), 'd:', round(res.x[8] * 1e6, 2), 'um')
    print('Total:', round((res.x[2] + res.x[5] + res.x[8]) * 1e6, 2), 'um')
    ns = array([res.x[0], res.x[3], res.x[6]])
    ks = array([res.x[1], res.x[4], res.x[7]])
    thikcs = array([res.x[2], res.x[5], res.x[8]])
else:
    print('Total:', round((res.x[2] + res.x[5]) * 1e6, 2), 'um')
    ns = array([res.x[0], res.x[3]])
    ks = array([res.x[1], res.x[4]])
    thikcs = array([res.x[2], res.x[5]])
print()
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()


H_teo = H_sim(ns, ks, thikcs, d_air, f_ref)


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
axs[0].set_xlim([f_ref[0], 1.2])
axs[0].set_ylim([0, 1])
axs[1].set_ylabel(r'$\phi \ (rad)$')
axs[1].set_xlim([f_ref[0], 1.2])
axs[1].set_ylim([- 4 * pi, pi])
xlabel(r'$f\ (THz)$')

E_sam_teo = irfft(H_teo * E_ref_w, n=t_sam.size)
figure(20)
plot(t_sam, E_sam, lw=1)
plot(t_sam, E_sam_teo, lw=1)
figure(21)
plot(t_sam, abs(E_sam - E_sam_teo), lw=1)


# fig, axs = subplots(2)
# fig.suptitle('Re/Im')
#
# axs[0].plot(f_ref, real(H_w), lw=1)
# axs[1].plot(f_ref, imag(H_w), lw=1)
#
# axs[0].plot(f_ref, real(H_teo), lw=1)
# axs[1].plot(f_ref, imag(H_teo), lw=1)
#
#
# axs[0].set_ylabel(r'$Re$')
# axs[0].xaxis.set_visible(False)
# axs[0].set_xlim([f_ref[0], f_ref[-1]])
# axs[1].set_ylabel(r'$Im$')
# axs[1].set_xlim([f_ref[0], f_ref[-1]])
# xlabel(r'$f\ (THz)$')


# figure(3)
# plot(t_ref, E_ref, lw=1, label='ref')
# plot(t_sam, E_sam, lw=1, label='sam')
# legend()
# figure(4)
# plot(f_ref, toDb(E_ref_w), lw=1, label='ref')
# plot(f_sam, toDb(E_sam_w), lw=1, label='sam')
# legend()


show()
