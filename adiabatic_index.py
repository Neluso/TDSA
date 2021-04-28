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


def get_H(M_mat):
    M_11 = M_mat[:, 0, 0]
    M_21 = M_mat[:, 1, 0]
    M_12 = M_mat[:, 0, 1]
    M_22 = M_mat[:, 1, 1]
    H_r = - M_21 / M_22
    H_t = (M_11 - M_12 * M_21 / M_22)
    # H_t = 1 / M_11
    return H_r, H_t


t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
delta_t_ref = mean(diff(t_ref))
enlargement = 1 * E_ref.size
ref_pulse_idx = centre_loc(E_ref)
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
D_adiab = 1e-4  # 100 um
n_1 = 1.4 + 1j * 0.03
thick_1 = 1e-3  # - D_adiab/2  # 1000 um
n_2 = 1.8 + 1j * 0.03
thick_2 = 1e-3  # - D_adiab/2  # 1000 um
phi_air = phase_factor(n_air, - thick_1 - thick_2, f_ref)
m = (n_2 - n_1) / D_adiab
b = n_1
# [1, 2, 10, 50, 100, 1000]:
for N_grid in [2, 50]:
    # N_grid = 10
    d = D_adiab / N_grid
    n_adiab = (arange(N_grid) + 0.5) * d
    n_adiab = f_n_1(n_1, n_2, D_adiab, n_adiab)
    n_adiab = append(n_1, n_adiab)
    n_adiab = append(n_adiab, n_2)
    d_adiab = d * ones(N_grid)
    d_adiab = append(thick_1 - D_adiab / 2, d_adiab)
    d_adiab = append(d_adiab, thick_2 - D_adiab / 2)

    M_i = list()
    N_i = list()

    for freq in f_ref:
        m_in = dot(P_i(n_adiab[0], d_adiab[0], freq), D_i_j(n_air, n_adiab[0]))
        n_mat_in = dot(P_i(n_1, thick_1, freq), D_i_j(n_air, n_1))
        m_i = m_in
        for i in range(N_grid):
            n_i = n_adiab[i]
            n_i_1 = n_adiab[i + 1]
            d_i = d_adiab[i + 1]
            m_j = dot(P_i(n_i, d_i, freq), D_i_j(n_i, n_i_1))
            m_i = dot(m_j, m_i)
        n_i = n_adiab[-2]
        n_i_1 = n_adiab[-1]
        d_i = d_adiab[-1]
        m_j = dot(P_i(n_i, d_i, freq), D_i_j(n_i, n_i_1))
        n_mat_j = dot(P_i(n_2, thick_2, freq), D_i_j(n_1, n_2))
        n_mat_i = dot(n_mat_j, n_mat_in)
        m_i = dot(m_j, m_i)
        n_mat_i = dot(D_i_j(n_2, n_air), n_mat_i)
        m_i = dot(D_i_j(n_i_1, n_air), m_i)
        M_i.append(m_i)
        N_i.append(n_mat_i)

    M_i = array(M_i)
    N_i = array(N_i)

    H_r, H_t = get_H(M_i)
    H_r *= phi_air
    H_t *= phi_air
    # H_t *= phi_air

    E_sim_r = irfft(H_r * E_ref_w)
    E_sim_t = irfft(H_t * E_ref_w)
    E_sim_r_w = H_r * E_ref_w
    E_sim_t_w = H_t * E_ref_w
    angle_r = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_r)
    angle_t = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_t)

    figure(1)
    title('R wave')
    plot(t_ref * 1e12, E_sim_r, label=str(N_grid), lw=1)
    legend()
    figure(2)
    title('T wave')
    plot(t_ref * 1e12, E_sim_t, label=str(N_grid), lw=1)
    legend()
    figure(3)
    title('R wave')
    plot(f_ref * 1e-12, angle_r, label=str(N_grid), lw=1)
    legend()
    figure(4)
    title('T wave')
    plot(f_ref * 1e-12, angle_t, label=str(N_grid), lw=1)
    legend()


H_r_2lay, H_t_2lay = get_H(N_i)
H_r_2lay *= phi_air
H_t_2lay *= phi_air
# H_t_2lay *= phi_air

E_sim_r_2lay = irfft(H_r_2lay * E_ref_w)
E_sim_t_2lay = irfft(H_t_2lay * E_ref_w)
E_sim_r_w_2lay = H_r_2lay * E_ref_w
E_sim_t_w_2lay = H_t_2lay * E_ref_w
angle_r_2lay = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_r_2lay)
angle_t_2lay = jepsen_unwrap(t_ref, E_ref, t_ref, E_sim_t_2lay)

figure(1)
# plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
title('R wave')
plot(t_ref * 1e12, E_sim_r_2lay, label='2_lay', lw=1)
legend()

figure(2)
title('T wave')
plot(t_ref * 1e12, E_ref, '--', label='ref', lw=1)
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
