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
    n_aux = real(n) - 1j * imag(n) * freq * 1e-12
    omg = 2 * pi * freq
    thick_aux = thick * cos(theta(n))
    phi = 2 * omg * thick_aux / c_0
    return exp(- 1j * n_aux * phi)


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


# t_ref, E_ref = read_1file('./ref_colorins.txt')
sample_num = 3
t_ref, E_ref = read_1file('./20210707_pegote/ref' + str(sample_num) + '.txt')
E_ref = - E_ref
t_sam, E_sam = read_1file('./20210707_pegote/sam' + str(sample_num) + '.txt')
# t_ref = arange(1, 5001) * 0.05e-12
# f_carrier = 0.3e12  # 300 Ghz
# E_ref = sin(2 * pi * t_ref * f_carrier)
# E_ref *= exp(- (t_ref - 50e-12)**2 / (2e-12)**2)
# t_ref *= 1e12

for sample_num in range(3):
    # t_ref, E_ref = read_1file('./20210707_capa_eff/sam' + str(sample_num + 1) + '.txt')
    # plot(t_ref, E_ref, label='capa_eff_'+str(sample_num+1), lw=1)
    t_ref, E_ref = read_1file('./20210707_pegote/sam' + str(sample_num + 1) + '.txt')
    plot(t_ref, E_ref, label='pegote_'+str(sample_num+1), lw=1)
legend()
show()
quit()

delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_ref[0] = 1
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_sam *= 1e-12
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_sam[0] = 1

D_adiab = 400e-6  # 100 um
n_1 = 1.45  # - 1j * 0.03  # blau
thick_1 = 800e-6  # - D_adiab/2  # 1000 um
n_2 = 1.55  # - 1j * 0.06  # groc
thick_2 = 600e-6  # - D_adiab/2  # 1000 um

phi_air = phase_factor(n_air, - thick_1 - thick_2 - D_adiab, f_ref)


# Building adiabatic samples
N_grid = 1


n_eff = (n_1 + n_2) / 2
H_teo_eff_adiab = -1  # cr_l_1_l(n_subs, n_2)
H_teo_eff_adiab = H_sim(f_ref, n_2, n_eff, thick_2, H_teo_eff_adiab)
H_teo_eff_adiab = H_sim(f_ref, n_eff, n_1, D_adiab, H_teo_eff_adiab)
H_teo_eff_adiab = H_sim(f_ref, n_1, n_air, thick_1, H_teo_eff_adiab)
H_teo_eff_adiab *= phi_air
E_sim_effad = irfft(H_teo_eff_adiab * E_ref_w)
H_teo_noad = -1  # cr_l_1_l(n_subs, n_2)
H_teo_noad = H_sim(f_ref, n_2, n_1, thick_2 + D_adiab/2, H_teo_noad)
H_teo_noad = H_sim(f_ref, n_1, n_air, thick_1 + D_adiab/2, H_teo_noad)
H_teo_noad *= phi_air
E_sim_noad = irfft(H_teo_noad * E_ref_w)
figure(1)
title('R wave')
t_ref *= 1e12
t_sam *= 1e12
# plot(t_ref, - E_ref, ',', label='ref', lw=1)
plot(t_sam, E_sam, '-.', label='sam', lw=1)
plot(t_ref, E_sim_effad, '--', label='eff', lw=1)
plot(t_ref, E_sim_noad, '--', label='noad', lw=1)
xlabel('t (ps)')
ylabel('Amplitud (u.a.)')
xlim([5, 45])
legend()
figure(2)
plot(t_sam, abs(irfft(E_sam_w / E_ref_w * wiener_filter(E_ref_w, beta=1e-5))), '-.', label='sam', lw=1)
plot(t_ref, abs(irfft(H_teo_eff_adiab)), '--', label='eff', lw=1)
plot(t_ref, abs(irfft(H_teo_noad)), '--', label='noad', lw=1)
legend()
show()
