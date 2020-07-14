from TDSA import *
import os
from scipy.optimize import differential_evolution, curve_fit

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


def espectral_guess(freq, amp, sigma, f0, sqr_pw):
    return amp * (10**(- ((freq - f0) / (2 * sigma))**2))


def smooth(M, span):
    M_aux = M
    p_max = M.size
    for p in range(p_max):
        if p - span < 0:
            M_aux[p] = sum(M[:2 * p + 1]) / (2 * p + 1)
        elif span < p - 1 < p_max - 1 - span:
            M_aux[p] = sum(M[p - span:p + span]) / (2 * span + 1)
        elif p + span > p_max - 1:
            M_aux[p] = sum(M[2 * p - p_max - 1:p_max - 1]) / (2 * p_max - 2 * p - 1)
    return M_aux


t0 = time_ns()
out_dir = './output/'
t_ref, E_ref = read_1file('./data/sim_resources/noise_ref.txt')  # t_ref in ps
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
wh = open(out_dir + 'resolution_limit.csv', 'a')

in_dir = './output/traces/'
dir_list = os.listdir(in_dir)

# if __name__ == '__main__':
for file_i in dir_list:
    t_sim, E_sim = read_1file(in_dir + file_i)  # t_ref in s
    f_sim, E_sim_w = fourier_analysis(t_sim, E_sim)  # f_ref in Hz

    d_mat, e_s_sim, e_inf_sim, tau_sim = file_i[:-4].split('_')
    d_mat = float(d_mat)
    e_s_sim = float(e_s_sim)
    e_inf_sim = float(e_inf_sim)
    tau_sim = float(tau_sim)

    f_min_idx, f_max_idx = f_min_max_idx(f_sim, 0, 20)

    # f_sim *= 1e-12   # f_ref in THz

    f_sim_cut = f_sim[f_min_idx:f_max_idx]
    E_sim_w_cut = E_sim_w[f_min_idx:f_max_idx]

    plot(f_sim, toDb_0(E_sim_w), lw=0.75)

    peak_list = find_peaks(abs(E_sim_w), height=0.9 * amax(abs(E_sim_w)))
    peak_idx = peak_list[0][0]

    p_bF = polyfit(f_sim_cut[peak_idx:], toDb_0(E_sim_w_cut)[peak_idx:], 1)
    m_bF = p_bF[0]
    b_bF = p_bF[1]  # m_bF: m  before fit, b_bF: b before fit

    for i in range(peak_idx + 1, f_sim_cut.size - 1):
        f_sim_4fit = f_sim_cut[i:]
        E_sim_w_4fit = toDb_0(E_sim_w_cut)[i:]
        p0 = polyfit(f_sim_4fit, E_sim_w_4fit, 1)
        delta_m = abs((m_bF - p0[0]) / m_bF)
        delta_b = abs((b_bF - p0[1]) / b_bF)

        # plot(i, delta_m, 'r.', i, delta_b, 'b.')

        if delta_m <= 0.001 and delta_b <= 0.001:
            plot(f_sim, p0[0] * f_sim + p0[1], lw=0.75)
            f_cutoff = f_sim[i]  # cut off freq in THz
            pDr = abs(E_sim_w[peak_idx]) - abs(E_sim_w[i])  # PEak Dynamic Range in dB
            break

        m_bF, b_bF = p0[0], p0[1]
    k_bounds = [
        (-100, 100e-6),  # d_air
        (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
        (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
        (0.99 * tau_sim, 1.01 * tau_sim),  # tau
        (0.1e-6, 1000e-6)  # d_mat
    ]
    # if d_mat < 1:
    #     k_bounds = [
    #         (-100, 100e-6),  # d_air
    #         (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
    #         (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
    #         (0.99 * tau_sim, 1.01 * tau_sim),  # tau
    #         (0.1e-6, 1e-6)  # d_mat
    #     ]
    # elif d_mat < 10:
    #     k_bounds = [
    #         (-100, 100e-6),  # d_air
    #         (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
    #         (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
    #         (0.99 * tau_sim, 1.01 * tau_sim),  # tau
    #         (1e-6, 10e-6)  # d_mat
    #     ]
    # elif d_mat < 100:
    #     k_bounds = [
    #         (-100, 100e-6),  # d_air
    #         (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
    #         (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
    #         (0.99 * tau_sim, 1.01 * tau_sim),  # tau
    #         (1e-6, 100e-6)  # d_mat
    #     ]
    # elif d_mat < 1000:
    #     k_bounds = [
    #         (-100, 100e-6),  # d_air
    #         (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
    #         (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
    #         (0.99 * tau_sim, 1.01 * tau_sim),  # tau
    #         (10e-6, 1000e-6)  # d_mat
    #     ]

    d_air_fit = list()
    e_s_fit = list()
    e_inf_fit = list()
    tau_fit = list()
    d_mat_fit = list()

    num_statistics = 10
    for i in range(num_statistics):
        print('Fitting', i + 1, 'of', num_statistics, 'for', d_mat, 'um')
        t1 = time_ns()
        res = differential_evolution(cost_function,
                                     k_bounds,
                                     args=(E_sim, E_ref_w, f_ref),
                                     popsize=90,
                                     maxiter=3000,
                                     # updating='deferred',
                                     # workers=-1,
                                     disp=False,  # step cost_function value
                                     polish=True
                                     )
        t2 = time_ns()
        # secs1 = (t2 - t1) * 1e-9
        # if secs1 < 3600:
        #     print('Fitting time (mm:ss):', strftime('%M:%S', gmtime(secs1)))
        # else:
        #     print('Fitting time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs1)))

        d_air_fit.append(res.x[0] * 1e6)
        e_s_fit.append(res.x[1])
        e_inf_fit.append(res.x[2])
        tau_fit.append(res.x[3])
        d_mat_fit.append(res.x[4] * 1e6)

    d_air_fit = array(d_air_fit)
    d_mat_fit = array(d_mat_fit)
    e_s_fit = array(e_s_fit)
    e_inf_fit = array(e_inf_fit)
    tau_fit = array(tau_fit)

    print('Saving simulation data for', d_mat, 'um')
    data = str(d_mat) + ',' + str(mean(d_mat_fit)) + ',' + str(std(d_mat_fit)) + ','
    data += str(e_s_sim) + ',' + str(mean(e_s_fit)) + ',' + str(std(e_s_fit)) + ','
    data += str(e_inf_sim) + ',' + str(mean(e_inf_fit)) + ',' + str(std(e_inf_fit)) + ','
    data += str(tau_sim) + ',' + str(mean(tau_fit)) + ',' + str(std(tau_fit)) + ','
    data += '0.0,' + str(mean(d_air_fit)) + ',' + str(std(d_air_fit)) + ','
    data += str(f_cutoff) + ',' + str(pDr) + '\n'

    wh.write(data)
    t3 = time_ns()
    secs0 = (t3 - t0) * 1e-9
    if secs0 < 3600:
        print('Time since start (mm:ss):', strftime('%M:%S', gmtime(secs0)))
    else:
        print('Time since start (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs0)))
    # quit()

wh.close()
print()
print('Finished')
t3 = time_ns()
secs0 = (t3 - t0) * 1e-9
if secs0 < 3600:
    print('Total processing time (mm:ss):', strftime('%M:%S', gmtime(secs0)))
else:
    print('Total processing time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs0)))