from TDSA import *
from scipy.optimize import differential_evolution


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
n_subs = 1.17e99  # substrate refractive index


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
    n_comp = n + 1j * k
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    return exp(2 * 1j * n_comp * omg * thick / c_0)


# def H_sim(n_l, k_l, thick_l, freq):
#     n_l_1 = n_aire
#     r_llp1 = 1
#     exp_beta = phase_factor(n_l, k_l, thick_l, freq)
#     r_lm1l = cr_l_1_l(n_l, n_l_1)
#     tt = ct2(n_l, n_l_1)
#     return r_lm1l + (tt * r_llp1 * exp_beta) / (1 + r_lm1l * r_llp1 * exp_beta)
#
#
# def cost_function(k, *args):
#     n1, k1, thick1 = k
#     H_meas, freqs = args
#     logs = log(abs(H_sim(n1, k1, thick1, freqs))) - log(abs(H_meas))
#     angles = angle(H_sim(n1, k1, thick1, freqs)) - angle(H_meas)
#     return sum(logs**2 + angles**2)
#     # return real(sum((H_sim(n1, k1, thick1, freqs) - H_meas)**2))

def H_sim(ns, ks, thicks, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, ns[-1]) * ones(freq.size)
    rlm1l = cr_l_1_l(ns[-1], ns[-2])
    tt = ct2(ns[-1], ns[-2])
    exp_phi = phase_factor(ns[-1], ks[-1], thicks[-1], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # inner layers, no constact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    for i in range(1, layers - 1):
        i = layers - i - 1
        rlm1l = cr_l_1_l(ns[i], ns[i - 1])
        tt = ct2(ns[i], ns[i - 1])
        exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
        H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # most outer layer, in contact with air
    rlm1l = cr_l_1_l(ns[0], n_air)
    tt = ct2(ns[0], n_air)
    exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


def cost_function(k, *args):
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
    H_meas, freqs = args
    logs = log(abs(H_sim(ns, ks, thicks, freqs))) - log(abs(H_meas))
    angles = unwrap(angle(H_sim(ns, ks, thicks, freqs))) - unwrap(angle(H_meas))
    return sum(logs**2 + angles**2)
    # return real(sum((H_sim(n1, k1, thick1, freqs) - H_meas)**2))

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

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))
f_min, f_max = f_min_max_idx(f_ref, 0.1, 1)
f_ref = f_ref[f_min:f_max]
E_ref_w = E_ref_w[f_min:f_max]
f_sam = f_sam[f_min:f_max]
E_sam_w = E_sam_w[f_min:f_max]


ref_pulse_idx = centre_loc(E_ref)
k_bounds = []

for i in range(layers):
    k_bounds.append((1, 5))  # n
    k_bounds.append((0, 1))  # k
    k_bounds.append((0, 1e-3))  # thickness

# TODO test different strategies
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam_w / E_ref_w, f_ref),
                             strategy='best1bin',
                             popsize=15
                             )

print(res)
plot(f_ref, abs(E_sam_w / E_ref_w), lw=1)
plot(f_ref, abs(H_sim(array([1.4, 1.5, 1.6]), 0.01 * ones(3), 33e-6 * ones(3), f_ref)), lw=1)
plot(f_ref, abs(H_sim(array([res.x[0], res.x[3], res.x[6]]), array([res.x[1], res.x[4], res.x[7]]), array([res.x[2], res.x[5], res.x[8]]), f_ref)), lw=1)
show()
