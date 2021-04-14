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


D_adiab = 1e-5  # 100 um
n_1 = 1.4 + 1j * 0.01
thick_1 = 1e-4  # 1 mm
n_2 = 1.6 + 1j * 0.01
thick_2 = 1e-4  # 1 mm
m = (n_2 - n_1) / D_adiab
b = n_1
N_grid = 2
d = D_adiab / N_grid
t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
M_i = list()

for freq in f_ref:
    m_in = dot(P_i(n_1, thick_1, freq), D_i_j(n_air, n_1))
    m_i = m_in
    for i in range(N_grid - 1):
        n_i = f_n_1(n_1, n_2, D_adiab, d * (i + 1))  # n_i
        n_i_1 = f_n_1(n_1, n_2, D_adiab, d * (i + 2))  # n_i+1
        phase_mat_i = P_i(n_i, d, freq)
        fres_mat_i = D_i_j(n_i, n_i_1)
        m_j = dot(phase_mat_i, fres_mat_i)
        m_i = dot(m_i, m_j)
    n_i = f_n_1(n_1, n_2, D_adiab, D_adiab - 1)  # n_i
    n_i_1 = n_2  # n_i+1
    phase_mat_i = P_i(n_i, thick_2, freq)
    fres_mat_i = D_i_j(n_i, n_i_1)
    m_out = dot(phase_mat_i, fres_mat_i)
    m_i = dot(m_i, m_out)
    M_i.append(m_i)
M_i = array(M_i)
M_11 = M_i[:, 0, 0]
M_21 = M_i[:, 1, 0]

phi_air = 1  # phase_factor(n_air, -(thick_1 + thick_2 + D_adiab), f_ref)
H_r = M_21 / M_11 * phi_air
H_t = 1 / M_11 * phi_air
E_sim = irfft(H_t * E_ref_w)
plot(t_ref, E_ref)
plot(t_ref, E_sim)
show()
