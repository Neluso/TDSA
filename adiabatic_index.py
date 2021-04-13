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


D_adiab = 1e-4  # 100 um
n_1 = 1.4
thick_1 = 1e-3  # 1 mm
n_2 = 1.6
thick_2 = 1e-3  # 1 mm
m = (n_2 - n_1) / D_adiab
b = n_1
t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
N_grid = 3
d = D_adiab / N_grid
M_in = dot(P_i(n_1, thick_1, f_ref), D_i_j(n_air, n_1))
# M_out = dot(P_i(n_2, thick_2, f_ref), D_i_j(n_2, n_air))
M_i = M_in
for i in range(N_grid - 1):
    n_i = f_n_1(n_1, n_2, D_adiab, d * (i + 1))  # n_i
    n_i_1 = f_n_1(n_1, n_2, D_adiab, d * (i + 2))  # n_i+1
    phase_mat_i = P_i(n_i, d, f_ref)
    fres_mat_i = D_i_j(n_i, n_i_1)
    m_i = dot(phase_mat_i, fres_mat_i)
    M_i = dot(M_i, m_i)
n_i = f_n_1(n_1, n_2, D_adiab, D_adiab - 1)  # n_i
n_i_1 = n_2  # n_i+1
phase_mat_i = P_i(n_i, thick_2, f_ref)
fres_mat_i = D_i_j(n_i, n_i_1)
M_out = dot(phase_mat_i, fres_mat_i)
M_i = dot(M_i, M_out)
