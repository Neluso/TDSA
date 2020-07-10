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


t0 = time_ns()
out_dir = './output/'

if __name__ == '__main__':
    t_ref, E_ref = read_1file('./data/sim_resources/transmision_ref.txt')  # t_ref in ps
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
    
    t_ref *= 1e-12  # t_ref in s
    f_ref *= 1e12   # f_ref in Hz
    
    
    # material data
    e_s_sim = 1.4**2
    e_inf_sim = 1.8**2
    tau_sim = 1e-12
    n_sim, k_sim = nk_from_eps(e_s_sim, e_inf_sim, tau_sim, f_ref)
    
    wh = open(out_dir + 'resolution_limit.csv', 'w')
    
    for d_mat in [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000]:
        
        print()
        print('Simulating for', d_mat, 'um')
        
        d_mat *= 1e-6
        d_air_sim = 0
        
        H_sim_teo = H_sim(f_ref, n_sim, k_sim, d_mat, d_air_sim)
        E_sim_w = H_sim_teo * E_ref_w
        E_sim = irfft(E_sim_w)
        
        name_trace = str(d_mat*1e6) + '_' + str(e_s_sim) + '_' + str(e_inf_sim) + '_' + str(tau_sim) + '.txt'
        print('Saving trace as', name_trace)
        tw = open(out_dir + name_trace, 'w')
        for i in range(E_sim.size):
            tw.write(str(t_ref[i]) + ',' + str(E_sim[i]) + '\n')
        tw.close()
        
        # k_bounds = [
        #     (0, 20e-6),  # d_air
        #     (1, 3 * e_s_sim + 1),  # e_s
        #     (1, 3 * e_inf_sim + 1),  # e_inf
        #     (0.0001 * tau_sim, 3 * tau_sim),  # tau
        #     (0, 3 * d_mat)  # d_mat
        # ]
    
        k_bounds = [
            (-100, 100e-6),       # d_air
            (1, 15),           # e_s
            (1, 15),           # e_inf
            (1e-15, 1e-9),     # tau
            (1e-6, 5 * d_mat)  # d_mat
        ]
        
        d_air_fit = list()
        e_s_fit = list()
        e_inf_fit = list()
        tau_fit = list()
        d_mat_fit = list()
        
        num_statistics = 10
        for i in range(num_statistics):
            print('Fitting', i + 1, 'of', num_statistics, 'for', d_mat * 1e6, 'um')
            t1 = time_ns()
            res = differential_evolution(cost_function,
                                         k_bounds,
                                         args=(E_sim, E_ref_w, f_ref),
                                         popsize=45,
                                         maxiter=3000,
                                         updating='deferred',
                                         workers=-1,
                                         disp=False,  # step cost_function value
                                         polish=True
                                         )
            t2 = time_ns()
            secs1 = (t2 - t1) * 1e-9
            if secs1 < 3600:
                print('Fitting time (mm:ss):', strftime('%M:%S', gmtime(secs1)))
            else:
                print('Fitting time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs1)))
            
            d_air_fit.append(res.x[0])
            e_s_fit.append(res.x[1])
            e_inf_fit.append(res.x[2])
            tau_fit.append(res.x[3])
            d_mat_fit.append(res.x[4])
    
        d_air_fit = array(d_air_fit)
        d_mat_fit = array(d_mat_fit)
        e_s_fit = array(e_s_fit)
        e_inf_fit = array(e_inf_fit)
        tau_sim = array(tau_sim)
        
        print('Saving simulation data for', d_mat * 1e6, 'um')
        data = str(d_mat) + ',' + str(mean(d_mat_fit)) + ',' + str(std(d_mat_fit)) + ','
        data += str(e_s_sim) + ',' + str(mean(e_s_fit)) + ',' + str(std(e_s_fit)) + ','
        data += str(e_inf_sim) + ',' + str(mean(e_inf_fit)) + ',' + str(std(e_inf_fit)) + ','
        data += str(tau_sim) + ',' + str(mean(tau_fit)) + ',' + str(std(tau_fit)) + ','
        data += str(d_air_sim) + ',' + str(mean(d_air_fit)) + ',' + str(std(d_air_fit)) + '\n'
        
        wh.write(data)
        t3 = time_ns()
        secs0 = (t3 - t0) * 1e-9
        if secs0 < 3600:
            print('Time since start (mm:ss):', strftime('%M:%S', gmtime(secs0)))
        else:
            print('Time since start (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs0)))
    
    wh.close()
    print()
    print('Finished')
    t3 = time_ns()
    secs0 = (t3 - t0) * 1e-9
    if secs0 < 3600:
        print('Total processing time (mm:ss):', strftime('%M:%S', gmtime(secs0)))
    else:
        print('Total processing time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs0)))
