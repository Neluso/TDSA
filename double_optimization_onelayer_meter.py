from TDSA import *
from scipy.optimize import differential_evolution
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air_cplx * sin(deg_in * pi / 180)
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index --- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index --- metal


# function definitions
def theta(n, theta_in):  # theta_in in degrees
    theta_in *= pi / 180
    return arcsin(n_air * theta_in / real(n))


def ct2(n_l, n_l_1, theta_in):
    n_l *= cos(theta(n_l, theta_in))
    n_l_1 *= cos(theta(n_l_1, theta_in))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1, theta_in):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1, theta_in))
    n_l *= cos(theta(n_l, theta_in))
    return - (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, k, thick, theta_in, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n, theta_in))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi) * exp(- k * phi)


def H_sim(n, k, thick, d_air, theta_in, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, n - 1j * k, theta_in) * ones(freq.size)
    rlm1l = cr_l_1_l(n - 1j * k, n_air_cplx, theta_in)
    tt = ct2(n - 1j * k, n_air_cplx, theta_in)
    exp_phi = phase_factor(n, k, thick, theta_in, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return phase_factor(n_air, 0, d_air, theta_in, freq) * H_i


def cost_function(params, *args):
    d_air, n, k, thick = params
    E_sam, E_ref_w, freqs, E_sam2, E_ref_w2, freqs2 = args
    H_teo = H_sim(n, k, thick, d_air, 30, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w)
    H_teo2 = H_sim(n, k, thick, d_air, 50, freqs2)
    E_sam_teo_w2 = E_ref_w2 * H_teo2
    E_sam_teo2 = irfft(E_sam_teo_w2)
    return sum((E_sam_teo - E_sam) ** 2 + (E_sam_teo2 - E_sam2) ** 2)


# Main script
# Boleto 176054
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat2_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_primer/ref metal primer_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_primer/sam metal primer_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref cork wcoat_avg_f.txt')

# 50º
# t_ref2, E_ref2 = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam2, E_sam2 = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/metal_w_coat/sam metal wcoat2_avg_f.txt')
# # t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# # t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/sam cork wcoat2_avg_f.txt')
# t_ref2, E_ref2 = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam2, E_sam2 = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/sam cork wcoat3_avg_f.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/sam cork wcoat 3_avg_f.txt')

# 50º
# t_ref2, E_ref2 = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017_50_deg/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam2, E_sam2 = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017_50_deg/cork_w_coat/sam cork wcoat3_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 3_avg_f.txt')
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/ref metal cork_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/sam cork wcoat 3_avg_f.txt')

# 50º
t_ref2, E_ref2 = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017_50_deg/cork_w_coat/ref metal cork wcoat_avg_f.txt')
t_sam2, E_sam2 = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017_50_deg/cork_w_coat/sam metal wcoat3_avg_f.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
ref_pulse_idx2 = centre_loc(E_ref2)
window = tukey(2 * ref_pulse_idx)
window2 = tukey(2 * ref_pulse_idx2)
window = zero_padding(window, 0, E_ref.size - window.size)
window2 = zero_padding(window2, 0, E_ref2.size - window2.size)
E_ref *= window
E_sam *= window
E_ref2 *= window2
E_sam2 *= window2
# E_ref = zero_padding(E_ref, 0, enlargement)
# t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
# E_sam = zero_padding(E_sam, 0, enlargement)
# t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
# E_sam_max = amax(abs(E_sam))


t_ref *= 1e-12
t_sam *= 1e-12
t_ref2 *= 1e-12
t_sam2 *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_ref2, E_ref_w2 = fourier_analysis(t_ref2, E_ref2)
f_sam2, E_sam_w2 = fourier_analysis(t_sam2, E_sam2)

k_bounds = []

k_bounds.append((0, 100e-6))     # air thickness
k_bounds.append((1, 5))          # n
k_bounds.append((0, 1))          # k
k_bounds.append((1e-6, 150e-6))  # thickness

n_mean = list()
k_mean = list()
thick_mean = list()
d_air_mean = list()
N_tests = 100
t1 = time_ns()
for i in trange(N_tests):
    res = differential_evolution(cost_function,
                                 k_bounds,
                                 args=(E_sam, E_ref_w, f_ref, E_sam2, E_ref_w2, f_ref2),
                                 disp=False,
                                 polish=True
                                 )
    # print(res)
    # print()
    # print('Avg layer')
    n_mean.append(res.x[1])
    k_mean.append(res.x[2])
    thick_mean.append(res.x[3])
    d_air_mean.append(res.x[0])
t2 = time_ns()

print('n =', round(mean(n_mean), 2), '+-', round(std(n_mean), 2))
print('k =', round(mean(k_mean), 3), '+-', round(std(k_mean), 3))
print('d =', round(mean(thick_mean) * 1e6, 2), '+-', round(std(thick_mean) * 1e6, 2), 'um')
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()

figure(50)
Ns_tests = arange(N_tests)
plot(Ns_tests, mean(n_mean) - n_mean, lw=1, label='n')
figure(51)
plot(Ns_tests, mean(k_mean) - k_mean, lw=1, label='k')
figure(52)
plot(Ns_tests, mean(thick_mean) - thick_mean, lw=1, label='t')
figure(53)
plot(Ns_tests, mean(d_air_mean) - d_air_mean, lw=1, label='d')
legend()
show()
quit()

H_w = E_sam_w / E_ref_w
H_teo = H_sim(mean(n_mean), mean(k_mean), mean(thick_mean), mean(d_air_mean), 30, f_ref)
H_teo2 = H_sim(mean(n_mean), mean(k_mean), mean(thick_mean), mean(d_air_mean), 50, f_ref2)
E_sam_teo = irfft(H_teo * E_ref_w)
E_sam_teo2 = irfft(H_teo2 * E_ref_w2)


t_ref *= 1e12
t_sam *= 1e12
t_ref2 *= 1e12
t_sam2 *= 1e12
f_ref *= 1e-12
f_sam *= 1e-12
f_ref2 *= 1e-12
f_sam2 *= 1e-12

figure(1)
plot(t_sam, E_sam, lw=1, label='30')
plot(t_sam2, E_sam2, lw=1, label='50')
plot(t_sam, E_sam_teo, lw=1, label='fit30')
plot(t_sam2, E_sam_teo2, lw=1, label='fit50')
legend()

fig, axs = subplots(2)
fig.suptitle('Abs/Phase')
axs[0].plot(f_ref, abs(H_w), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_w)), lw=1)
axs[0].plot(f_ref, abs(H_teo), lw=1)
axs[1].plot(f_ref, unwrap(angle(H_teo)), lw=1)
axs[0].set_ylabel(r'$\rho$')
# axs[0].xaxis.set_visible(False)
axs[1].set_ylabel(r'$\phi \ (rad)$')
xlabel(r'$f\ (THz)$')
# axs[0].set_xlim([f_ref[0], 1.2])
# axs[0].set_ylim([0, 1])
# axs[1].set_xlim([f_ref[0], 1.2])
# axs[1].set_ylim([- 4 * pi, pi])

show()
