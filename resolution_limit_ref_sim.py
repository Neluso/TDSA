from TDSA import *
from scipy import stats
from scipy.optimize import curve_fit


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


def abs_spec(freq, s, sigma=1, mu=0):
    x = (freq - mu) / sigma
    norm_fact = sqrt(2 * pi * log(10))
    y = (log10(x)**2) / (2 * s**2)
    exp = 10**(-y) / (x * s * norm_fact)
    # exp = 10 ** (-y) / (x * s)
    return exp


def sim_refs():
    out_dir = './output/refs/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        
    for trash_file in os.listdir(out_dir):
        os.remove(out_dir + trash_file)
    
    # num_points = 1251
    num_points = 2501
    # num_points = 10001
    freqs = arange(num_points) * 1e10  # Hz
    freqs *= 1000 / (num_points - 1)
    # print(freqs)
    # quit()
    freqs *= 1e-12  # THz
    
    # print(mean(diff(freqs)))
    # quit()
    
    times = arange(2 * (num_points - 1))
    times = times / (mean(diff(freqs)) * 2 * (num_points - 1))  # ps
    
    E_sim_w_abs = stats.lognorm.pdf(freqs, 0.5, scale=0.4)
    # E_sim_w_abs = stats.lognorm.pdf(freqs, 1, scale=0.4, loc=0.15)
    E_sim_w_abs2 = stats.lognorm.pdf(freqs, 1, scale=0.4, loc=0.158)
    E_sim_w_abs3 = stats.lognorm.pdf(freqs, 0.2, scale=0.4, loc=-0.08)
    E_sim_w_abs4 = stats.lognorm.pdf(freqs, 1.5, scale=1, loc=0.2)
    
    E_sim_w_abs /= max(E_sim_w_abs)
    E_sim_w_abs2 /= max(E_sim_w_abs2)
    E_sim_w_abs3 /= max(E_sim_w_abs3)
    E_sim_w_abs4 /= max(E_sim_w_abs4)
    
    # E_sim_w_abs = fromDb(30 + toDb(E_sim_w_abs))
    E_sim_w_arg = phase_factor(n_air, 0, 1.6e-3, freqs*1e12)
    # E_sim_w_arg2 = phase_factor(n_air, 0, 1.6e-3, freqs * 1e12)
    # E_sim_w_arg3 = phase_factor(n_air, 0, 1.6e-3, freqs * 1e12)
    # E_sim_w_arg4 = phase_factor(n_air, 0, 1.6e-3, freqs * 1e12)
    
    # plot(freqs, E_sim_w_arg)
    E_sim_ref = irfft(E_sim_w_abs * E_sim_w_arg)  # * 100
    E_sim_ref2 = irfft(E_sim_w_abs2 * E_sim_w_arg)  # * 100
    E_sim_ref3 = irfft(E_sim_w_abs3 * E_sim_w_arg)  # * 100
    E_sim_ref4 = irfft(E_sim_w_abs4 * E_sim_w_arg)  # * 100
    # E_sim_ref /= max(abs(E_sim_ref))
    # E_sim_ref2 /= max(abs(E_sim_ref2))
    # E_sim_ref3 /= max(abs(E_sim_ref3))
    
    # for ns_floor in [-90, -70, -60, -40, -30, -20, -10]:
    # for ns_floor in [-90, -60, -30, -10]:
    for ns_floor in [-60]:  # -90, -60, -40]:
        num_traces = 10
        trace_statitics = zeros(E_sim_ref.shape)
        trace_statitics2 = zeros(E_sim_ref2.shape)
        trace_statitics3 = zeros(E_sim_ref3.shape)
        trace_statitics4 = zeros(E_sim_ref4.shape)
        for j in range(num_traces):
            # E_sim_ref += fromDb(ns_floor) * random.normal(0, 0.01, E_sim_ref.size)
            trace_statitics += E_sim_ref + fromDb(ns_floor) * random.normal(0, 0.02, E_sim_ref.size)
            trace_statitics2 += E_sim_ref2 + fromDb(ns_floor) * random.normal(0, 0.02, E_sim_ref2.size)
            trace_statitics3 += E_sim_ref3 + fromDb(ns_floor) * random.normal(0, 0.02, E_sim_ref3.size)
            trace_statitics4 += E_sim_ref4 + fromDb(ns_floor) * random.normal(0, 0.02, E_sim_ref4.size)
        if ns_floor == -60:
            figure()
            plot(freqs, toDb_0(rfft(trace_statitics / num_traces)), lw=1)  # , label='old')
            # plot(freqs, toDb_0(rfft(trace_statitics2 / num_traces)), lw=1, label='new')
            # plot(freqs, toDb_0(rfft(trace_statitics3 / num_traces)), lw=1, label='newer')
            # plot(freqs, toDb_0(rfft(trace_statitics4 / num_traces)), lw=1, label='newest')
            # legend()
            xlabel(r'$f\ (THz)$')
            ylabel('dB')
            savefig('./output/ref_spectra_newest.png')
            close()
            figure()
            plot(times, 100 * trace_statitics / num_traces, lw=1)  # , label='old')
            # plot(times, 100 * trace_statitics2 / num_traces, lw=1, label='new')
            # plot(times, 100 * trace_statitics3 / num_traces, lw=1, label='newer')
            # plot(times, 100 * trace_statitics4 / num_traces, lw=1, label='newest')
            # legend()
            xlabel(r'$t\ (ps)$')
            savefig('./output/ref_traces.png')
            close()
        write_data(times, 100 * trace_statitics / num_traces, str(ns_floor) + '_ref', out_dir)  # THz


sim_refs()
show()
