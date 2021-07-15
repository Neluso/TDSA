from TDSA import *


deg_in = 45  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 2.6  # n_air_cplx


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def ct2(n_l, n_l_1):
    n_l_aux = n_l * cos(theta(n_l))
    n_l_1_aux = n_l_1 * cos(theta(n_l_1))
    return 4 * n_l_aux * n_l_1_aux / (n_l_aux + n_l_1_aux)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_aux = n_l * cos(theta(n_l))
    n_l_1_aux = n_l_1 * cos(theta(n_l_1))
    return (n_l_1_aux - n_l_aux) / (n_l_1_aux + n_l_aux)


def phase_factor(n, thick, freq):  # theta in radians
    n = real(n) - 1j * imag(n) * freq * 1e-12
    omg = 2 * pi * freq
    thick_aux = thick * cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi)


def H_sim(freq, n_l, n_l_1, thick, H_subs):
    # print(n_l, n_l_1)
    n_l = real(n_l) - 1j * imag(n_l) * freq * 1e-12
    n_l_1 = real(n_l_1) - 1j * imag(n_l_1) * freq * 1e-12
    H_i = H_subs  # cr_l_1_l(n_subs, n - 1j * k)
    rlm1l = cr_l_1_l(n_l, n_l_1)
    tt = ct2(n_l, n_l_1)
    exp_phi = phase_factor(n_l, thick, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


def f_n_1(n_1, n_2, D_adiab, d):
    m = (n_2 - n_1) / D_adiab
    b = n_1
    return m * d + b  # linar transition
    # return m * d**2 / D_adiab + b  # quadratic transition
    # return m * d**3 / D_adiab**2 + b  # cubic transition


def path_p_num_plot(file_path, sample_num):
    t_ref, E_ref = read_1file(file_path + 'ref' + str(sample_num) + '.txt')
    E_ref = - E_ref
    t_sam_ad, E_sam_ad = read_1file(file_path + 'sam' + str(sample_num) + '.txt')
    t_sam_nad, E_sam_nad = read_1file(file_path + 'sam' + str(sample_num) + '.txt')
    delta_t_ref = mean(diff(t_ref))
    enlargement = 0 * E_ref.size
    ref_pulse_idx = centre_loc(E_ref)
    E_ref = zero_padding(E_ref, 0, enlargement)
    t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    t_ref *= 1e-12
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_ref[0] = 1
    E_sam_ad = zero_padding(E_sam_ad, 0, enlargement)
    t_sam_ad = concatenate((t_sam_ad, t_sam_ad[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    t_sam_ad *= 1e-12
    f_sam_ad, E_sam_ad_w = fourier_analysis(t_sam_ad, E_sam_ad)
    f_sam_ad[0] = 1
    E_sam_nad = zero_padding(E_sam_nad, 0, enlargement)
    t_sam_nad = concatenate((t_sam_nad, t_sam_nad[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    t_sam_nad *= 1e-12
    f_sam_nad, E_sam_nad_w = fourier_analysis(t_sam_nad, E_sam_nad)
    f_sam_nad[0] = 1
    H_ad_w = E_sam_ad_w / E_ref_w
    # H_ad_w = E_ref_w / E_ref_w
    H_nad_w = E_sam_nad_w / E_ref_w
    # plot(f_ref * 1e-12, toDb(E_sam_w))
    # plot(f_ref * 1e-12, toDb_0(E_ref_w), lw=1)
    # xlim([0, 4])
    # ylim([-100, 5])
    # xlabel(r'$f\ (THz)$')
    # ylabel(r'$(dB)$')
    # plot(f_ref * 1e-12, jepsen_unwrap(t_ref, E_ref, t_sam, E_sam))
    # plot(f_ref * 1e-12, unwrap(angle(E_sam_w)))
    plot(t_sam_ad * 1e12, E_sam_ad / amax(abs(E_sam_ad)), lw=1, label='sam')
    # plot(t_ref * 1e12, E_ref / amax(abs(E_ref)), lw=1)
    deconv_ad = irfft(H_ad_w * wiener_filter(E_ref_w, beta=0.01))**2
    plot(t_sam_ad * 1e12, roll(deconv_ad, centroid_E2(t_ref, E_ref)) / amax(abs(deconv_ad)), lw=1, label=r'$dcv^2$')
    # xlim([0, 50])
    # plot(t_sam_nad * 1e12, E_sam_nad / amax(abs(E_sam_nad)), lw=1)
    # deconv_nad = abs(irfft(H_nad_w * wiener_filter(E_ref_w, beta=1e-2)))
    # plot(t_sam_nad * 1e12, deconv_nad / amax(abs(deconv_nad)), lw=1)
    # plot(t_ref * 1e12, - E_ref, lw=1)
    # plot(t_sam_ad * 1e12, E_sam_ad, lw=1)
    # plot(t_sam_nad * 1e12, E_sam_nad, lw=1)
    # plot(f_ref * 1e-12, toDb(E_ref_w * wiener_filter(E_ref_w, beta=1e-5)))
    # plot(E_sam)


time_ref = 180
t_ref1, E_ref1 = read_1file('./' + str(time_ref) + 's_1.txt')
t_ref2, E_ref2 = read_1file('./' + str(time_ref) + 's_2.txt')
# for time_ref in [10, 20, 30, 50, 90, 120, 180]:
#     t_ref1, E_ref1 = read_1file('./' + str(time_ref) + 's_1.txt')
#     t_ref2, E_ref2 = read_1file('./' + str(time_ref) + 's_2.txt')
#     ctr1_idx = centroid_E2(t_ref1, E_ref1)
#     ctr2_idx = centroid_E2(t_ref2, E_ref2)
#     figure(3)
#     plot(time_ref, t_ref1[ctr1_idx], 'b.')
#     plot(time_ref, t_ref2[ctr2_idx], 'r.')
#     figure(1)
#     plot(t_ref1, E_ref1, label=str(time_ref), lw=1)
#     figure(2)
#     plot(t_ref2, E_ref2, label=str(time_ref), lw=1)
# figure(1)
# xlim([20, 30])
# legend()
# figure(2)
# xlim([20, 30])
# legend()
# show()
# quit()
E_ref1 *= -1
E_ref2 *= -1
delta_t_ref = mean(diff(t_ref1))
enlargement = 0 * E_ref1.size
E_ref1 = zero_padding(E_ref1, 0, enlargement)
t_ref1 = concatenate((t_ref1, t_ref1[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref1 *= 1e-12
E_ref2 = zero_padding(E_ref2, 0, enlargement)
t_ref2 = concatenate((t_ref2, t_ref2[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref2 *= 1e-12
f_ref1, E_ref1_w = fourier_analysis(t_ref1, E_ref1)
f_ref1[0] = 1
f_ref2, E_ref2_w = fourier_analysis(t_ref2, E_ref2)
f_ref2[0] = 1

D_adiab = 50e-6  # 100 um
n_1 = 1.45  # - 1j * 0.03  # blau
thick_1 = 1000e-6  # - D_adiab/2  # 1000 um
n_2 = 1.55  # - 1j * 0.06  # groc
thick_2 = 1000e-6  # - D_adiab/2  # 1000 um
phi_air = phase_factor(n_air, - thick_1 - thick_2, f_ref1)
N_grid = 1000

for D_adiab in [1000e-6]:  # , 1e-4, 1.5e-4]:  # , 1e-3]:  # 1000 és suficient per a simular adiabàtic
    # "d" adaptable, N_grid fixe
    # d = D_adiab / N_grid
    # "d" fixe, N_grid adaptable
    d = 1e-7  # 100 nm
    N_grid = int(D_adiab / d)
    m = (n_2 - n_1) / D_adiab
    b = n_1
    n_adiab = [n_1]
    for i in range(N_grid):
        n_adiab.append(f_n_1(n_1, n_2, D_adiab, d * i))
    n_adiab.append(n_2)
    d_adiab = d * ones(N_grid)
    d_adiab = append(thick_1 - D_adiab / 2, d_adiab)
    d_adiab = append(d_adiab, thick_2 - D_adiab / 2)

    n_adiab = append(n_air_cplx, n_adiab)
    d_adiab = append(0, d_adiab)

    H_teo_adiab = -1  # cr_l_1_l(n_subs, n_adiab[-1])

    n_layers = d_adiab.size
    for i in range(1, n_layers):
        j = n_layers - i
        H_teo_adiab = H_sim(f_ref1, n_adiab[j], n_adiab[j - 1], d_adiab[j], H_teo_adiab)

    H_teo_adiab *= phi_air
    E_sim_ad = irfft(H_teo_adiab * E_ref1_w)
    noise_fig = E_ref1_w / E_ref2_w
    figure(time_ref)
    plot(t_ref1 * 1e12, abs(irfft(H_teo_adiab * noise_fig)), lw=0.5, label='noisy')
    plot(t_ref1 * 1e12, abs(irfft(H_teo_adiab * noise_fig * wiener_filter(E_ref2_w, beta=1e-2))), lw=1, label='filt')
    plot(t_ref1 * 1e12, abs(irfft(H_teo_adiab)), lw=0.8, label='pure')
    # xlim([10, 25])
    # ylim([0, 0.002])
    legend()


# plot(irfft(H_teo_adiab))
# plot(irfft(E_ref1_w / E_ref2_w), lw=1)
# plot(irfft(E_ref2_w / E_ref1_w), lw=1)

show()
