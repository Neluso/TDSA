from TDSA import *


deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = n_air_cplx


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
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi)


def H_sim(freq, n_l, n_l_1, thick, H_subs):
    H_i = H_subs  # cr_l_1_l(n_subs, n - 1j * k)
    rlm1l = cr_l_1_l(n_l, n_l_1)
    tt = ct2(n_l, n_l_1)
    exp_phi = phase_factor(n_l, thick, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


def f_n_1(n_1, n_2, D_adiab, d):
    m = (n_2 - n_1) / D_adiab
    b = n_1
    return m * d + b


t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
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
D_adiab = 5e-3  # 100 um
n_1 = 1.4 - 1j * 0.03
thick_1 = 1e-4  # - D_adiab/2  # 1000 um
n_2 = 2.6 - 1j * 0.03
thick_2 = 1e-4  # - D_adiab/2  # 1000 um
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
N_grid = 10000
# for D_adiab in [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]:  # 1000 és suficient per a simular adiabàtic
# for D_adiab in [1e-6, 10**-5.5, 1e-5, 10**-4.5, 1e-4, 10**-3.5, 1e-3]:
for D_adiab in [1e-6, 1e-5, 1e-4]:  # , 1e-4, 1.5e-4]:  # , 1e-3]:  # 1000 és suficient per a simular adiabàtic
    # # "d" adaptable, N_grid fixe
    # d = D_adiab / N_grid
    # # "d" fixe, N_grid adaptable
    d = 1e-10
    N_grid = int(D_adiab / d)
    n_adiab = (arange(N_grid) + 0.5) * d
    n_adiab = f_n_1(n_1, n_2, D_adiab, n_adiab)
    n_adiab = append(n_1, n_adiab)
    n_adiab = append(n_adiab, n_2)
    d_adiab = d * ones(N_grid)
    d_adiab = append(thick_1 - D_adiab / 2, d_adiab)
    d_adiab = append(d_adiab, thick_2 - D_adiab / 2)

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
    plot(t_ref * 1e12, E_sim_ad, label=legend_Text, lw=1)
    figure(2)
    plot(f_ref * 1e-12, abs(H_teo_adiab), label=legend_Text, lw=1)
    figure(3)
    plot(t_ref * 1e12, irfft(H_teo_adiab), label=legend_Text, lw=1)
    figure(4)
    plot(t_ref * 1e12, real(irfft(H_teo_adiab - H_teo_no_adiab)), label=legend_Text, lw=1)
    print(D_adiab * 1e6, 'um ->', sum(abs(irfft(H_teo_adiab - H_teo_no_adiab))) * 1e3)
    figure(5)
    plot(D_adiab * 1e6, sum(abs(irfft(H_teo_adiab - H_teo_no_adiab))), 'b.')


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
plot(t_ref * 1e12, E_sim_nad, '--', label='2_lay', lw=1)
# plot(t_ref * 1e12, E_sim_effad, '-.', label='2_lay_eff', lw=1)
legend()
figure(2)
plot(f_ref * 1e-12, abs(H_teo_no_adiab), '--', label='2_lay', lw=1)
# plot(f_ref * 1e-12, abs(H_teo_eff_adiab), '-.', label='2_lay_eff', lw=1)
legend()
figure(3)
plot(t_ref * 1e12, irfft(H_teo_no_adiab), '--', label='2_lay', lw=1)
# plot(t_ref * 1e12, abs(irfft(H_teo_eff_adiab)), '-.', label='2_lay_eff', lw=1)
legend()
figure(4)
legend()
figure(5)
xlabel('D_adiab')
ylabel('sum diff H')
show()
