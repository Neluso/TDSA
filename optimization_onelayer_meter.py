from TDSA import *
from scipy.optimize import differential_evolution


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air_cplx * sin(deg_in * pi / 180)
n_subs = 1.17e99 - 0.0 * 1j  # substrate refractive index


# function definitions
def theta(n):
    return arcsin(snell_sin / n)


def ct2(n_l, n_l_1):
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    thetal_1 = theta(n_l_1)
    thetal = theta(n_l)
    return (n_l_1 * cos(thetal_1) - n_l * cos(thetal)) / (n_l_1 * cos(thetal_1) + n_l * cos(thetal))


def phase_factor(n, k, thick, freq):  # theta in radians
    n_comp = n - 1j * k
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    # return exp(- 2 * 1j * n_comp * omg * thick / c_0)
    return exp(1j * n * phi) * exp(- k * phi)


def H_sim(n, k, thick, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, n - 1j * k) * ones(freq.size)
    rlm1l = cr_l_1_l(n - 1j * k, n_air_cplx)
    tt = ct2(n - 1j * k, n_air_cplx)
    exp_phi = phase_factor(n, k, thick, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # # inner layers, no constact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    # for i in range(1, layers - 1):
    #     i = layers - i - 1
    #     rlm1l = cr_l_1_l(ns[i], ns[i - 1])
    #     tt = ct2(ns[i], ns[i - 1])
    #     exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
    #     H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    #
    # # most outer layer, in contact with air
    # rlm1l = cr_l_1_l(ns[0], n_air_cplx)
    # tt = ct2(ns[0], n_air_cplx)
    # exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    # H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


def cost_function(params, *args):
    n, k, thick = params
    H_meas, freqs = args
    logs = log(abs(H_sim(n, k, thick, freqs))) - log(abs(H_meas))
    angles = abs(unwrap(angle(H_sim(n, k, thick, freqs)))) - abs(unwrap(angle(H_meas)))
    return sum(logs + angles)


# Main script
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')

delta_t_ref = mean(diff(t_ref))
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
delta_f_ref = mean(diff(f_ref))
f_min, f_max = f_min_max_idx(f_ref, 0.2, 0.9)
f_ref = f_ref[f_min:f_max]
E_ref_w = E_ref_w[f_min:f_max]
f_sam = f_sam[f_min:f_max]
E_sam_w = E_sam_w[f_min:f_max]


ref_pulse_idx = centre_loc(E_ref)
k_bounds = []


k_bounds.append((1, 5))  # n
k_bounds.append((0, 1))  # k
k_bounds.append((0, 1e-3))  # thickness

H_w = E_sam_w / E_ref_w
# H_w /= amax(H_w)

res = differential_evolution(cost_function,
                             k_bounds,
                             args=(H_w, f_ref),
                             strategy='best1bin',
                             popsize=150,
                             maxiter=2500,
                             )
# print(res)
# print()
print('Avg layer')
print('n =', round(res.x[0], 2))
print('k =', round(res.x[1], 3))
print('d =', round(res.x[2] * 1e6, 2), 'um')
print()

H_teo = H_sim(res.x[0],
              res.x[1],
              res.x[2],
              f_ref)


f_ref *= 1e-12

fig, axs = subplots(2)
fig.suptitle('Abs/Phase')

axs[0].plot(f_ref, abs(H_w), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_w)), lw=1)

axs[0].plot(f_ref, abs(H_teo), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_teo)), lw=1)


fig.legend(labels=['sample', 'fit'])
axs[0].set_ylabel(r'$\rho$')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], f_ref[-1]])
axs[1].set_ylabel(r'$\phi$')
axs[1].set_xlim([f_ref[0], f_ref[-1]])
xlabel(r'$f\ (THz)$')


show()
