from TDSA import *
from resolution_limit_ref_sim import sim_refs


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
    phi = 2 * omg * thick / c_0
    return exp(- 1j * phi) * exp(- k * phi)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air):
    
    H_i = cr_l_1_l(n_subs, n_i - 1j * k_i)
    
    rlm1l = cr_l_1_l(n_i - 1j * k_i, n_o - 1j * k_o)
    tt = ct2(n_i - 1j * k_i, n_o - 1j * k_o)
    exp_phi = phase_factor(n_i, k_i, thick_i, freq)
    
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_o - 1j * k_o, n_air_cplx)
    tt = ct2(n_o - 1j * k_o, n_air_cplx)
    exp_phi = phase_factor(n_o, k_o, thick_o, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air, e_s_i, e_inf_i, tau_i, thick_i, e_s_o, e_inf_o, tau_o, thick_o = params
    E_sam, E_ref_w, freqs = args
    n_i, k_i = nk_from_eps(e_s_i, e_inf_i, tau_i, freqs)  # debye model
    n_o, k_o = nk_from_eps(e_s_o, e_inf_o, tau_o, freqs)  # debye model
    H_teo = H_sim(freqs, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air)
    E_teo = irfft(H_teo * E_ref_w)
    return sum((E_sam - E_teo)**2)


def sim_traces():
    t0 = time_ns()
    out_dir = './output/traces/'
    in_dir = './output/refs/'
    ref_list = os.listdir(in_dir)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    for trash_file in os.listdir(out_dir):
        os.remove(out_dir + trash_file)
    
    for ref_file in ref_list:
        t_ref, E_ref = read_1file(in_dir + ref_file)  # t_ref in ps
        f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
        
        t_ref *= 1e-12  # t_ref in s
        f_ref *= 1e12   # f_ref in Hz
        
        ns_level = ref_file.split('_')[0]
        
        # material data
        e_s_sim_i = 1.4**2
        e_inf_sim_i = 1.8**2
        tau_sim_i = 5e-14
        n_sim_i, k_sim_i = nk_from_eps(e_s_sim_i, e_inf_sim_i, tau_sim_i, f_ref)
        e_s_sim_o = 1.6 ** 2
        e_inf_sim_o = 1.8 ** 2
        tau_sim_o = 5e-14
        n_sim_o, k_sim_o = nk_from_eps(e_s_sim_o, e_inf_sim_o, tau_sim_o, f_ref)
        
        f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0, 1)
        f_ref *= 1e-12  # THz
        figure()
        plot(f_ref, n_sim_i, f_ref, n_sim_o)
        xlabel(r'$f\ (THz)$')
        # xlim([f_ref[f_min_idx], f_ref[f_max_idx]])
        # ylim([0.9 * n_sim[f_min_idx], 1.1 * n_sim[f_max_idx]])
        savefig('./output/n_sim.png')
        close()
        
        figure()
        plot(f_ref, k_sim_i, f_ref, k_sim_o)
        xlabel(r'$f\ (THz)$')
        # xlim([f_ref[f_min_idx], f_ref[f_max_idx]])
        # ylim([k_sim[f_min_idx], k_sim[f_max_idx]])
        savefig('./output/k_sim.png')
        close()
        # plot(t_ref*1e12, - E_ref, label='ref')
        
        # plot(f_ref*1e12, unwrap(angle(E_ref_w)), label='ref')
        # for d_mat in [0.1, 0.2, 0.3, 0.4, 0.5]:
        # for d_mat in [1, 2, 3, 4, 5]:
        # for d_mat in [10, 20, 30, 40, 50]
        # for d_mat in [100, 200, 300, 400, 500]:
        # for d_mat in [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500]:
        # for d_mat in [0.1, 0.5, 1, 5, 10, 50, 100, 500]:
        # for d_mat in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100]:
        # for d_mat in [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1.1e3]:
        # for d_mat in [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3]:
        
        f_ref *= 1e12  # Hz
        for d_mat in [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500]:
            
            print()
            print('Simulating for', d_mat, 'um')
            name_trace = str(d_mat).zfill(6) + '_' + str(e_s_sim_i) + '_' + str(e_inf_sim_i) + '_'
            name_trace = name_trace + str(tau_sim_i) + '_' + str(e_s_sim_o) + '_' + str(e_inf_sim_o) + '_'
            name_trace = name_trace + str(tau_sim_o) + '_' + ns_level + '.txt'

            d_mat *= 1e-6  # um
            
            H_sim_teo = H_sim(f_ref, n_sim_i, k_sim_i, d_mat, n_sim_o, k_sim_o, d_mat, 0)  # - d_mat)
            # plot(f_ref, unwrap(angle(H_sim_teo)))
            # show()
            # quit()
            E_sim_w = H_sim_teo * E_ref_w
            E_sim = irfft(E_sim_w)
            # plot(t_ref*1e12, E_sim, label=round(d_mat*1e6, 1))
            # plot(f_ref, unwrap(angle(E_sim_w)), label=round(d_mat * 1e6, 1))
            
            print('Saving trace as', name_trace)
            tw = open(out_dir + name_trace, 'w')
            for i in range(E_sim.size):
                tw.write(str(t_ref[i]*1e12) + ',' + str(E_sim[i]) + '\n')
            tw.close()
        # legend()
        # show()


print('Simulating refs to "measure" samples')
print()
sim_refs()
print('Simulating traces - "measuring" samples')
sim_traces()
print()
print('Simulating refs to "measure" references')
sim_refs()
print()
print('----------------------------')
print('Done')
print('----------------------------')
