from TDSA import *
from scipy.optimize import differential_evolution, NonlinearConstraint
from scipy.signal.windows import tukey

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air_cplx * sin(deg_in * pi / 180)
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index --- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index --- metal


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


def f_cons(params):
    d_air, n, k, thick_1 = params
    # thick = thick_1
    thick = 1 / thick_1
    return abs(n * thick * cos(theta(abs(n))))


def H_sim(n, k, thick, d_air, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, n - 1j * k) * ones(freq.size)
    rlm1l = cr_l_1_l(n - 1j * k, n_air_cplx)
    tt = ct2(n - 1j * k, n_air_cplx)
    exp_phi = phase_factor(n, k, thick, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return phase_factor(n_air, 0, d_air, freq) * H_i


def cost_function(params, *args):
    d_air, n, k, thick_1 = params
    E_sam, E_ref_w, freqs = args
    thick = 1 / thick_1
    # thick = thick_1
    H_teo = H_sim(n, k, thick, d_air, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    return sum((E_sam_teo - E_sam) ** 2)


# Main script
# Boleto 176054
# @ 30º
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_primer/ref metal primer_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_primer/sam metal primer_avg_f.txt')
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref cork wcoat_avg_f.txt')
# @ 50º
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018_50_deg/cork_w_coat/sam cork wcoat2_avg_f.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/sam cork wcoat 3_avg_f.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
window = tukey(2 * ref_pulse_idx)
window = zero_padding(window, 0, E_ref.size - window.size)
E_ref *= window
E_sam *= window
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))
E_sam_max_idx = where(abs(E_sam) == E_sam_max)[0][0]

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

H_w = E_sam_w / E_ref_w
H_w_filt = H_w * wiener_filter(E_ref_w, beta=fromDb(noise_floor(f_ref, E_ref_w, 20)))
irf_filt = irfft(H_w_filt, t_ref.size)
irf_filt_max = amax(irf_filt)
irf_filt_max_idx = where(irf_filt == irf_filt_max)[0][0]
irf = roll(irf_filt / amax(abs(irf_filt)), E_sam_max_idx - irf_filt_max_idx)  # align deconvolution to t_ref
# plot(t_ref, irf)
# show()
# quit()
irf_dno = SWT_denoising(irf, 3, 0.1)  # apply de-noising
irf_peaks = signal.find_peaks(abs(irf), 0.3 * amax(abs(irf)))  # peak to be, at least, 30% of maximum
# irf_peaks = signal.find_peaks(abs(irf_dno), 0.3 * amax(abs(irf_dno)))  # peak to be, at least, 30% of maximum
print(irf_peaks)
# times = diff(t_ref[irf_peaks[0]])
# times = diff(t_ref[array([1892, 2085])])
# times = diff(t_ref[array([1894, 2062])])  # cork 1
# times = diff(t_ref[array([1903, 2079])])  # cork 2
times = diff(t_ref[array([1863, 2037])])  # cork 3


k_bounds = []
# testing
k_bounds.append((0, 100e-6))     # air thickness
k_bounds.append((1, 5))          # n
k_bounds.append((0, 1))          # k
# k_bounds.append((1e-6, 150e-6))  # thickness
# k_bounds.append((1e4, 1e6))  # 1/thickness
k_bounds.append((0, 1e20))  # 1/thickness for testing

cons_error = 0.1

# cons = NonlinearConstraint(f_cons, 0.8 * times[0]*c_0/2, 1.2 * times[0]*c_0/2#,
#                            jac='2-point',
#                            hess='2-point'
#                            )
cons = NonlinearConstraint(f_cons, (1 - cons_error) * times[0]*c_0/2, (1 + cons_error) * times[0]*c_0/2)
# cons = LinearConstraint(array([0, 1, 0, - times[0]*c_0/2]), -0.1, 0.1)

res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             disp=True,
                             constraints=cons,
                             polish=True
                             )
res.x[3] = 1 / res.x[3]
print(res)
# S = zeros((4, 4))
# S = dot(res.jac[0],res.jac[1])
# S = res.jac[0]
# S = dot(S.T, S)
# S = pinv(S)

print()
# print('Avg layer')
# print('n =', round(res.x[1], 2), '+-', round(sqrt(S[1, 1]), 2))
# print('k =', round(res.x[2], 3), '+-', round(sqrt(S[2, 2]), 3))
# print('d =', round(res.x[3] * 1e6, 2), '+-', round(sqrt(1/S[3, 3]) * 1e6, 2), 'um')
print('Avg layer')
print('n =', round(res.x[1], 2))
print('k =', round(res.x[2], 3))
print('d =', round(res.x[3] * 1e6, 2), 'um')
print()

H_teo = H_sim(res.x[1], res.x[2], res.x[3], res.x[0], f_ref)
E_sam_teo = irfft(H_teo * E_ref_w, n=t_sam.size)
print('Fit goodnes:', sum(abs(E_sam - E_sam_teo)))

t_sam *= 1e12
t_ref *= 1e12
f_ref *= 1e-12
f_sam *= 1e-12

figure(1)
# plot(t_ref, E_ref, lw=1, label='ref')
plot(t_sam, E_sam, lw=1, label='sam')
plot(t_sam, E_sam_teo, lw=1, label='fit')
xlim([t_ref[0], t_ref[-1]])
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
