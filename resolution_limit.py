from TDSA import *
from scipy.optimize import differential_evolution


deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal


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


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi_n = 2 * omg * thick / c_0
    phi_k = 2 * omg * thick / c_0
    return exp(- 1j * phi_n) * exp(- k * phi_k)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, n, k, thick, d_air):
    
    H_i = cr_l_1_l(n_subs, n - 1j * k)
    rlm1l = cr_l_1_l(n - 1j * k, n_air_cplx)
    tt = ct2(n - 1j * k, n_air_cplx)
    exp_phi = phase_factor(n, k, thick, freq)
    
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air, e_s, e_inf, tau, thick = params
    E_sam, E_ref_w, freqs = args
    n, k = nk_from_eps(e_s, e_inf, tau, freqs)  # debye model
    H_teo = H_sim(freqs, n, k, thick, d_air)
    E_teo = irfft(H_teo * E_ref_w)
    return sum((E_sam - E_teo)**2)


t_ref, E_ref = read_1file('./data/sim_resources/transmision_ref.txt')  # t_ref in ps
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz

t_ref *= 1e-12  # t_ref in s
f_ref *= 1e12  # f_ref in Hz


# material data
e_s_sim = 2**2
e_inf_sim = 2.2**2
tau_sim = 1e-14
n_sim, k_sim = nk_from_eps(e_s_sim, e_inf_sim, tau_sim, f_ref)


# d_mat = 100e-6
for d_mat in tqdm([1e-6, 10e-6, 50e-6, 100e-6, 150e-6, 500e-6, 1e-3]):
    H_sim_teo = H_sim(f_ref, n_sim, k_sim, d_mat, 0)
    E_sim_w = H_sim_teo * E_ref_w
    E_sim = irfft(E_sim_w)
    
    # k_bounds = [
    #     (0, 20e-6),  # d_air
    #     (1, 3 * e_s_sim + 1),  # e_s
    #     (1, 3 * e_inf_sim + 1),  # e_inf
    #     (0.0001 * tau_sim, 3 * tau_sim),  # tau
    #     (0, 3 * d_mat)  # d_mat
    # ]

    k_bounds = [
        (0, 100e-6),  # d_air
        (1, 15),  # e_s
        (1, 15),  # e_inf
        (1e-15, 1e-9),  # tau
        (1e-6, 5 * d_mat)  # d_mat
    ]
    
    res = differential_evolution(cost_function,
                                 k_bounds,
                                 args=(E_sim, E_ref_w, f_ref),
                                 popsize=45,
                                 maxiter=3000,
                                 disp=False,  # step cost_function value
                                 polish=True
                                 )
    
    # print(res)
    d_air_fit = res.x[0]
    e_s_fit = res.x[1]
    e_inf_fit = res.x[2]
    tau_fit = res.x[3]
    d_mat_fit = res.x[4]
    
    figure(50)
    plot(d_mat*1e6, d_mat_fit*1e6, 'ro')
    title('Fitted vs actual thickness')
    xlabel('d_mat (um)')
    ylabel('d_mat_fit (um)')
    figure(51)
    plot(d_mat*1e6, 100 * abs(d_mat-d_mat_fit) / d_mat, 'ro')
    xlabel('d_mat (um)')
    ylabel('delta_d')
    title('Error vs actual thickness')

show()
quit()

#
# print('Fitted values')
# print('d_mat =', d_mat_fit*1e6, '- err', 100 * abs(d_mat-d_mat_fit) / d_mat, '%')
# print('e_s =', e_s_fit, '- err', 100 * abs(e_s_fit-e_s_sim) / e_s_sim, '%')
# print('e_inf =', e_s_fit, '- err', 100 * abs(e_inf_fit-e_inf_sim) / e_inf_sim, '%')
# print('tau =', tau_fit, '- err', 100 * abs(tau_fit-tau_sim) / tau_sim, '%')
#
# n_fit, k_fit = nk_from_eps(e_s_fit, e_inf_fit, tau_fit, f_ref)
# H_fit = H_sim(f_ref, n_fit, k_fit, d_mat_fit, d_air_fit)
# E_fit_w = H_fit * E_ref_w
# E_fit = irfft(E_fit_w)
#
# figure(1)
# # plot(t_ref, E_ref, label='ref')
# plot(t_ref, E_sim, label='sim')
# plot(t_ref, E_fit, label='fit')
# legend()
#
# # figure(2)
# # plot(f_ref, toDb(E_ref_w), label='ref')
# # plot(f_ref, toDb(E_sam_w), label='sam')
# # legend()
#
# figure(3)
# plot(f_ref, n_sim, label='sim')
# plot(f_ref, n_fit, label='fit')
# title('n')
# legend()
#
# figure(4)
# plot(f_ref, k_sim, label='sim')
# plot(f_ref, k_fit, label='fit')
# title('k')
# legend()
#
# show()
