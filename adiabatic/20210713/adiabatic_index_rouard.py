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
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1))
    n_l *= cos(theta(n_l))
    return (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, thick, freq):  # theta in radians
    n = real(n) - 1j * imag(n) * freq * 1e-12
    omg = 2 * pi * freq
    thick *= cos(theta(n))
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
    t_ref_aux, E_ref_aux = read_1file(file_path + 'ref' + str(sample_num + 1) + '.txt')
    E_ref -= E_ref_aux
    t_ref_aux, E_ref_aux = read_1file(file_path + 'ref' + str(sample_num + 2) + '.txt')
    E_ref -= E_ref_aux
    t_sam_ad, E_sam_ad = read_1file(file_path + 'sam' + str(sample_num) + '.txt')
    t_ref_aux, E_sam_ad_aux = read_1file(file_path + 'sam' + str(sample_num + 1) + '.txt')
    E_sam_ad += E_sam_ad_aux
    t_ref_aux, E_sam_ad_aux = read_1file(file_path + 'sam' + str(sample_num + 2) + '.txt')
    E_sam_ad += E_sam_ad_aux
    E_ref /= 3
    E_sam_ad /= 3
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
    H_ad_w = abs(E_sam_ad_w) / abs(E_ref_w)
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
    # plot(t_sam_ad * 1e12, E_sam_ad / amax(abs(E_sam_ad)), lw=1)  # , label='sam')
    # plot(t_sam_ad * 1e12, E_sam_ad, lw=1, label='sam')
    # plot(t_ref * 1e12, E_ref / amax(abs(E_ref)), lw=1)
    deconv_ad = irfft(H_ad_w * wiener_filter(E_ref_w, beta=0.01))
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


# sample_num = 3
# file_path = './20210712_pegote/'
# file_path = ./20210712_pegote/'
# file_path = './20210712_3_capas_adiab/'
# path_p_num_plot('./20210712_capa_eff/', 1)
# path_p_num_plot('./20210712_3_capas_adiab/', 1)
# figure(1)
# # path_p_num_plot('./20210712_vidre/20210712_adiab/', 1)
# path_p_num_plot('./20210713_adiab/3_capas/', 1)
# # path_p_num_plot('./20210713_adiab/3_capas/', 2)
# # path_p_num_plot('./20210713_adiab/3_capas/', 3)
# path_p_num_plot('./20210715_adiab/3_capas/', 1)
# # path_p_num_plot('./20210715_adiab/3_capas/', 2)
# # path_p_num_plot('./20210715_adiab/3_capas/', 3)
# title('3 capas')
# legend()
# figure(2)
# path_p_num_plot('./20210713_adiab/capa_eff/', 1)
# # path_p_num_plot('./20210713_adiab/capa_eff/', 2)
# # path_p_num_plot('./20210713_adiab/capa_eff/', 3)
# path_p_num_plot('./20210715_adiab/capa_eff/', 1)
# # path_p_num_plot('./20210715_adiab/capa_eff/', 2)
# # path_p_num_plot('./20210715_adiab/capa_eff/', 3)
# path_p_num_plot('./20210715_adiab/capa_eff_2/', 1)
# path_p_num_plot('./20210715_adiab/capa_eff_3/', 1)
# # hlines(-50, 0, 4, linestyles='dashed', colors='r', label='-50 dB')
# # hlines(-60, 0, 4, linestyles='dashed', colors='k', label='-60 dB')
# # hlines(-70, 0, 4, linestyles='dashed', colors='b', label='-70 dB')
# # hlines(-60, 0, 4, linestyles='dashed', label='-60 dB')
# vlines(13, -2, 2, linestyles='dashed', colors='y', label='amarillo')
# vlines(20, -2, 2, linestyles='dashed', colors='g', label='verde')
# vlines(29, -2, 2, linestyles='dashed', colors='b', label='azul')
# vlines(34.7, -2, 2, linestyles='dashed', colors='k', label='substrato')
# title('capa eff')
# xlabel(r'$t\ (ps)$')
# legend()
figure(3)
path_p_num_plot('./20210715_spray_adiab/', 1)
# path_p_num_plot('./20210715_spray_adiab/', 2)
# path_p_num_plot('./20210715_spray_adiab/', 3)
# vlines(22.5, -2, 2, linestyles='dashed', colors='r', label='rojo')
# vlines(23, -2, 2, linestyles='dashed', colors='g', label='verde')
# vlines(23.6, -2, 2, linestyles='dashed', colors='b', label='plastico')
# vlines(25.2, -2, 2, linestyles='dashed', colors='k', label='substrato')
# vlines(27.9, -2, 2, linestyles='dashed', colors='b', label='subs-plas2')

legend()


# figure(2)
# path_p_num_plot('./20210712_vidre/20210712_adiab/', 2)


show()
