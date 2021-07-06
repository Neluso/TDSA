from TDSA import *


deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
n_subs = 2.6  # n_air_cplx


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


# t_ref, E_ref = read_1file('./ref_colorins.txt')
t_ref, E_ref = read_1file('./ref.txt')
# t_ref = arange(1, 5001) * 0.05e-12
# f_carrier = 0.3e12  # 300 Ghz
# E_ref = sin(2 * pi * t_ref * f_carrier)
# E_ref *= exp(- (t_ref - 50e-12)**2 / (2e-12)**2)
# t_ref *= 1e12


delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_ref[0] = 1
# print(f_ref * 1e-12)
# quit()
D_adiab = 50e-6  # 100 um
n_1 = 1.45 - 1j * 0.03  # blau
thick_1 = 5e-5  # - D_adiab/2  # 1000 um
n_2 = 1.55 - 1j * 0.06  # groc
thick_2 = 5e-5  # - D_adiab/2  # 1000 um


freq_aux, n_b, n_b_std, alpha_b, alpha_b_std = read_from_1file('./blava.txt')
n_b = interp(f_ref, freq_aux, n_b, left=n_b[0], right=n_b[-1])
alpha_b = np.interp(f_ref, freq_aux, alpha_b, left=alpha_b[0], right=alpha_b[-1])
k_b = 1e-10 * c_0 * alpha_b / (4 * pi * f_ref)


freq_aux, n_g, n_g_std, alpha_g, alpha_g_std = read_from_1file('./groga.txt')
n_g = interp(f_ref, freq_aux, n_g, left=n_g[0], right=n_g[-1])
alpha_g = np.interp(f_ref, freq_aux, alpha_g, left=alpha_g[0], right=alpha_g[-1])
k_g = 1e-10 * c_0 * alpha_g / (4 * pi * f_ref)
# phi_air = phase_factor(n_air, - thick_1 - thick_2, f_ref)
phi_air = 1
m = (n_2 - n_1) / D_adiab
b = n_1

# Two layer sample
H_teo_no_adiab = cr_l_1_l(n_subs, n_2)
H_teo_no_adiab = H_sim(f_ref, n_2, n_1, thick_2, H_teo_no_adiab)
H_teo_no_adiab = H_sim(f_ref, n_1, n_air, thick_1, H_teo_no_adiab)
H_teo_no_adiab *= phi_air
E_sim_nad = irfft(H_teo_no_adiab * E_ref_w)

# [1, 2, 10, 50, 100, 1000]:
# Building adiabatic samples
N_grid = 4
# for D_adiab in [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]:  # 1000 és suficient per a simular adiabàtic
# for D_adiab in [1e-6, 10**-5.5, 1e-5, 10**-4.5, 1e-4, 10**-3.5, 1e-3]:
# for D_adiab in [1e-5, 10**-4.9, 10**-4.8, 10**-4.7, 10**-4.6, 10**-4.5, 10**-4.4, 10**-4.3, 10**-4.2, 10**-4.1, 1e-4, 10**-3.9, 10**-3.8, 10**-3.7, 10**-3.6, 10**-3.5, 10**-3.4, 10**-3.3, 10**-3.2, 10**-3.1, 1e-3]:
for D_adiab in [200e-6]:   # , 1e-4, 5e-4]:  # , 1e-4, 1.5e-4]:  # , 1e-3]:  # 1000 és suficient per a simular adiabàtic
    # "d" adaptable, N_grid fixe
    d = D_adiab / N_grid
    # "d" fixe, N_grid adaptable
    # d = 1e-7  # 100 nm
    # N_grid = int(D_adiab / d)
    m = (n_2 - n_1) / D_adiab
    b = n_1
    n_adiab = [n_1]
    for i in range(N_grid):
        n_adiab.append(f_n_1(n_1, n_2, D_adiab, d * i))
    n_adiab.append(n_2)
    print(len(n_adiab))


    print(d * 1e6, 'um')
    print('Grid:', N_grid)
    # n_adiab = (arange(N_grid) + 0.5) * d
    # n_adiab = f_n_1(n_1, n_2, D_adiab, n_adiab)
    # n_adiab = append(n_1, n_adiab)
    # n_adiab = append(n_adiab, n_2)
    # n_adiab *= 0.01 * (2 * random.rand(N_grid + 2) - 1)
    d_adiab = d * ones(N_grid)
    d_adiab = append(thick_1 - D_adiab / 2, d_adiab)
    d_adiab = append(d_adiab, thick_2 - D_adiab / 2)
    # d_adiab = append(thick_1, d_adiab)  # - D_adiab / 2, d_adiab)
    # d_adiab = append(d_adiab, thick_2)  # - D_adiab / 2)

    n_adiab = append(n_air_cplx, n_adiab)
    d_adiab = append(0, d_adiab)

    H_teo_adiab = cr_l_1_l(n_subs, n_adiab[-1])

    n_layers = d_adiab.size
    for i in range(1, n_layers):
        j = n_layers - i
        H_teo_adiab = H_sim(f_ref, n_adiab[j], n_adiab[j - 1], d_adiab[j], H_teo_adiab)

    H_teo_adiab *= phi_air
    E_sim_ad = irfft(H_teo_adiab * E_ref_w)

    legend_Text = str(D_adiab * 1e6) + ' ' + 'um'
    figure(1)
    # print(t_ref.size)
    # quit()
    plot(t_ref * 1e12, E_sim_ad, label=legend_Text, lw=1)
    # figure(2)
    # plot(f_ref * 1e-12, toDb(H_teo_adiab * E_ref_w), label=legend_Text, lw=1)
    figure(3)
    plot(t_ref * 1e12, irfft(H_teo_adiab), label=legend_Text, lw=1)
    # figure(4)
    # plot(t_ref * 1e12, real(irfft(H_teo_adiab - H_teo_no_adiab)), label=legend_Text, lw=1)
    # print(D_adiab * 1e6, 'um ->', sum(abs(irfft(H_teo_adiab - H_teo_no_adiab))) * 1e3)
    # figure(5)
    # plot(D_adiab * 1e6, sum(abs(irfft(H_teo_adiab - H_teo_no_adiab))), 'b.')


n_eff = (n_1 + n_2) / 2
H_teo_eff_adiab = cr_l_1_l(n_subs, n_2)
H_teo_eff_adiab = H_sim(f_ref, n_2, n_eff, thick_2 - D_adiab / 2, H_teo_eff_adiab)
H_teo_eff_adiab = H_sim(f_ref, n_eff, n_1, D_adiab, H_teo_eff_adiab)
H_teo_eff_adiab = H_sim(f_ref, n_1, n_air, thick_1 - D_adiab / 2, H_teo_eff_adiab)
H_teo_eff_adiab *= phi_air
E_sim_effad = irfft(H_teo_eff_adiab * E_ref_w)
figure(1)
# plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
title('R wave')
t_ref *= 1e12
plot(t_ref, E_sim_nad, '--', label='2_lay', lw=1)
plot(t_ref * 1e12, E_sim_effad, '-.', label='2_lay_eff', lw=1)
xlabel('t (ps)')
ylabel('Amplitud (u.a.)')
# xlim([t_ref[0], t_ref[-1]])
xlim([5, 45])
legend()
# figure(2)
# plot(f_ref * 1e-12, toDb(H_teo_eff_adiab * E_ref_w), '--', label='2_lay', lw=1)
# plot(f_ref * 1e-12, abs(H_teo_eff_adiab), '-.', label='2_lay_eff', lw=1)
legend()
figure(3)
plot(t_ref, irfft(H_teo_no_adiab), '--', label='2_lay', lw=1)
# xlim([1, 4])
# ylim([-0.005, 0.001])
# # plot(t_ref * 1e12, abs(irfft(H_teo_eff_adiab)), '-.', label='2_lay_eff', lw=1)
legend()
# figure(4)
# # legend()
# # figure(5)
# # xlabel('D_adiab')
# # ylabel('sum diff H')
show()
