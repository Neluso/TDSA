# Criteri de signes ñ = n + 1j * k


from TDSA import *


def t_i_j(n_i, n_j):
    return (2 * n_i) / (n_i + n_j)


def r_i_j(n_i, n_j):
    return (n_i - n_j) / (n_i + n_j)


def phase_factor(n_i, thick_i, freq):  # theta in radians
    omg = 2 * pi * freq
    phi = omg * thick_i / c_0
    return exp(1j * n_i * phi)


def P_i(n_i, thick_i, freq):
    phi_fac_1 = phase_factor(n_i, thick_i, freq)
    phi_fac_2 = phase_factor(n_i, - thick_i, freq)
    return array(((phi_fac_1, 0), (0, phi_fac_2)))


def D_i_j(n_i, n_j):
    t_12 = t_i_j(n_i, n_j)
    inv_t_12 = 1 / t_12
    return inv_t_12 * array(((1, r_i_j(n_i, n_j)), (r_i_j(n_i, n_j), 1)))


def f_n_1(n_1, n_2, D_adiab, d):
    m = (n_2 - n_1) / D_adiab
    b = n_1
    return m * d + b


t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# [1, 2, 10, 50, 100, 1000]:
for N_grid in [10, 50]:
    D_adiab = 1e-5  # 10 um
    n_1 = 1.4 + 1j * 0.03
    thick_1 = 1e-4 - D_adiab/2  # 100 um
    n_2 = 1.6 + 1j * 0.03
    thick_2 = 1e-4 - D_adiab/2  # 100 um
    m = (n_2 - n_1) / D_adiab
    b = n_1
    # N_grid = 10
    d = D_adiab / N_grid

    M_i = list()
    N_i = list()

    for freq in f_ref:
        m_in = dot(P_i(n_1, thick_1, freq), D_i_j(n_air, n_1))
        m_i = m_in
        for i in range(N_grid):
            n_i = f_n_1(n_1, n_2, D_adiab, d * (i + 1))  # n_i
            n_i_1 = f_n_1(n_1, n_2, D_adiab, d * (i + 2))  # n_i+1
            phase_mat_i = P_i(n_i, d, freq)
            fres_mat_i = D_i_j(n_i, n_i_1)
            m_j = dot(phase_mat_i, fres_mat_i)
            m_i = dot(m_i, m_j)
        n_i = n_2  # f_n_1(n_1, n_2, D_adiab, D_adiab - d)  # n_i
        n_i_1 = n_air  # n_i+1
        phase_mat_i = P_i(n_i, thick_2, freq)
        fres_mat_i = D_i_j(n_i, n_i_1)
        m_out = dot(phase_mat_i, fres_mat_i)
        m_i = dot(m_i, m_out)
        n_i = dot(m_out, m_in)
        M_i.append(m_i)
        N_i.append(n_i)

    M_i = array(M_i)
    M_11 = M_i[:, 0, 0]
    M_21 = M_i[:, 1, 0]
    M_12 = M_i[:, 0, 1]
    M_22 = M_i[:, 1, 1]

    N_i = array(N_i)
    N_11 = N_i[:, 0, 0]
    N_21 = N_i[:, 1, 0]
    N_12 = N_i[:, 0, 1]
    N_22 = N_i[:, 1, 1]

    phi_air = phase_factor(n_air, - (thick_1 + thick_2 + D_adiab), f_ref)
    H_r = - M_21 / M_22 * phi_air
    H_t = (M_11 - M_12*M_21 / M_22) * phi_air
    # H_r_2lay = - M_21 / M_22 * phi_air
    # H_t_2lay = (M_11 - M_12 * M_21 / M_22) * phi_air

    E_sim_r = irfft(H_r * E_ref_w)
    E_sim_t = irfft(H_t * E_ref_w)
    E_sim_r_w = H_r * E_ref_w
    E_sim_t_w = H_t * E_ref_w
    angle_r = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_r)
    angle_t = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_t)

    figure(1)
    # plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
    title('R wave')
    plot(t_ref * 1e12, E_sim_r, label=str(N_grid), lw=1)
    legend()
    figure(2)
    title('T wave')
    plot(t_ref * 1e12, E_sim_t, label=str(N_grid), lw=1)
    legend()
    figure(3)
    # plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
    title('R wave')
    plot(f_ref * 1e-12, angle_r, label=str(N_grid), lw=1)
    legend()
    figure(4)
    title('T wave')
    plot(f_ref * 1e-12, angle_t, label=str(N_grid), lw=1)
    legend()
    # plot(t_ref * 1e12, E_ref / amax(abs(E_ref)))
    # plot(t_ref * 1e12, E_sim / amax(abs(E_sim)))

    # figure(2)
    # plot(f_ref * 1e-12, abs(H_r))
    # plot(f_ref * 1e-12, abs(H_t))

H_r_2lay = - M_21 / M_22 * phi_air
H_t_2lay = (M_11 - M_12 * M_21 / M_22) * phi_air

E_sim_r_2lay = irfft(H_r * E_ref_w)
E_sim_t_2lay = irfft(H_t * E_ref_w)
E_sim_r_w_2lay = H_r * E_ref_w
E_sim_t_w_2lay = H_t * E_ref_w
angle_r_2lay = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_r_2lay)
angle_t_2lay = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_t_2lay)

figure(1)
# plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
title('R wave')
plot(t_ref * 1e12, E_sim_r_2lay, label='2_lay', lw=1)
legend()
figure(2)
title('T wave')
plot(t_ref * 1e12, E_sim_t_2lay, label='2_lay', lw=1)
legend()
figure(3)
# plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
title('R wave')
plot(f_ref * 1e-12, angle_r_2lay, label='2_lay', lw=1)
legend()
figure(4)
title('T wave')
plot(f_ref * 1e-12, angle_t_2lay, label='2_lay', lw=1)
legend()

show()
