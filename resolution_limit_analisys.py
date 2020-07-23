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
in_refs = './output/refs/'

dir_list = os.listdir(in_dir)

if __name__ == '__main__':
    for file_i in dir_list:
        t_sim, E_sim = read_1file(in_dir + file_i)  # t_ref in ps
        f_sim, E_sim_w = fourier_analysis(t_sim, E_sim)  # f_ref in THz
        
        ns_level = file_i.split('_')[-1]
        ns_level_file = ns_level.replace('.', '_ref.')
        
        t_ref, E_ref = read_1file(in_refs + ns_level_file)  # t_ref in ps
        f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
    
        d_mat_str, e_s_sim_str, e_inf_sim_str, tau_sim_str, pDr_ref = file_i[:-4].split('_')
        
        d_mat = float(d_mat_str)
        e_s_sim = float(e_s_sim_str)
        e_inf_sim = float(e_inf_sim_str)
        tau_sim = float(tau_sim_str)
        
        f_min_idx, f_max_idx = f_min_max_idx(f_sim*1e12, 0.35, 1.5)
        p1 = polyfit(f_sim[f_min_idx:f_max_idx], toDb_0(E_sim_w[f_min_idx:f_max_idx]), 1)
        m, b = p1[0], p1[1]
        f_cutoff = (float(ns_level.split('.')[0]) - b) / m  # THz
        # plot(f_sim, toDb_0(E_sim_w))
        # plot(f_sim, m * f_sim + b)
        
        # k_bounds = [  # 1% uncertainty in optical paramaters
        #         (-1e-12, 1e-12),  # d_air
        #         (0.99 * e_s_sim, 1.01 * e_s_sim),  # e_s
        #         (0.99 * e_inf_sim, 1.01 * e_inf_sim),  # e_inf
        #         (0.99 * tau_sim, 1.01 * tau_sim),  # tau
        #         (0, 2e-3)  # d_mat
        # ]
        # k_bounds = [  # 5% uncertainty in optical paramaters
        #     (-1e-12, 1e-12),  # d_air
        #     (0.95 * e_s_sim, 1.05 * e_s_sim),  # e_s
        #     (0.95 * e_inf_sim, 1.05 * e_inf_sim),  # e_inf
        #     (0.95 * tau_sim, 1.05 * tau_sim),  # tau
        #     (0, 2e-3)  # d_mat
        # ]
        k_bounds = [  # 10% uncertainty in optical paramaters
            (-1e-12, 1e-12),  # d_air
            (0.9 * e_s_sim, 1.1 * e_s_sim),  # e_s
            (0.9 * e_inf_sim, 1.1 * e_inf_sim),  # e_inf
            (0.9 * tau_sim, 1.1 * tau_sim),  # tau
            (0, 2e-3)  # d_mat
        ]
        # k_bounds = [  # 15% uncertainty in optical paramaters
        #     (-1e-12, 1e-12),  # d_air
        #     (0.85 * e_s_sim, 1.15 * e_s_sim),  # e_s
        #     (0.85 * e_inf_sim, 1.15 * e_inf_sim),  # e_inf
        #     (0.85 * tau_sim, 1.15 * tau_sim),  # tau
        #     (0, 2e-3)  # d_mat
        # ]
        # k_bounds = [  # 50% uncertainty in optical paramaters
        #     (-1e-12, 1e-12),  # d_air
        #     (0.5 * e_s_sim, 1.5 * e_s_sim),  # e_s
        #     (0.5 * e_inf_sim, 1.5 * e_inf_sim),  # e_inf
        #     (0.5 * tau_sim, 1.5 * tau_sim),  # tau
        #     (0, 1e-3)  # d_mat
        # ]
        # k_bounds = [  # 100% uncertainty in optical paramaters
        #     (-1e-12, 1e-12),  # d_air
        #     (0.05 * e_s_sim, 2 * e_s_sim),  # e_s
        #     (0.05 * e_inf_sim, 2 * e_inf_sim),  # e_inf
        #     (0.05 * tau_sim, 2 * tau_sim),  # tau
        #     (0, 2e-3)  # d_mat
        # ]
        # k_bounds = [  # >100% uncertainty in optical paramaters
        #     (-1e-12, 1e-12),  # d_air
        #     (1, 3 * e_s_sim),  # e_s
        #     (1, 3 * e_inf_sim),  # e_inf
        #     (0, 1),  # tau
        #     (0, 2e-3)  # d_mat
        # ]
        d_air_fit = list()
        e_s_fit = list()
        e_inf_fit = list()
        tau_fit = list()
        d_mat_fit = list()
        
        f_ref *= 1e12  # Hz
        f_sim *= 1e12  # Hz
        num_statistics = 10
        for i in range(num_statistics):
            print('Fitting', i + 1, 'of', num_statistics, 'for', d_mat, 'um at', ns_level.split('.')[0], 'dB')
            t1 = time_ns()
            res = differential_evolution(cost_function,
                                         k_bounds,
                                         args=(E_sim, E_ref_w, f_ref),
                                         popsize=90,
                                         maxiter=3000,
                                         updating='deferred',
                                         workers=-1,
                                         disp=True,  # step cost_function value
                                         polish=True
                                         )
            t2 = time_ns()
            print(res)
            n_sim, k_sim = nk_from_eps(res.x[1], res.x[2], res.x[3], f_sim)
            H_fit = H_sim(f_sim, n_sim, k_sim, res.x[4], res.x[0])
            # plot(t_sim, E_sim)
            # plot(t_sim, irfft(H_fit * E_ref_w))
            # show()
            # quit()
            secs1 = (t2 - t1) * 1e-9
            if secs1 < 3600:
                print('Fitting time (mm:ss):', strftime('%M:%S', gmtime(secs1)))
            else:
                print('Fitting time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs1)))
            t3 = time_ns()
            secs0 = (t3 - t0) * 1e-9
            if secs0 < 3600:
                print('Time since start (mm:ss):', strftime('%M:%S', gmtime(secs0)))
            else:
                print('Time since start (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs0)))
    
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
        data += str(f_cutoff) + ',' + str(pDr_ref) + '\n'
    
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
