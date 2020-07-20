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
    for trash_file in os.listdir(out_dir):
        os.remove(out_dir + trash_file)
    
    num_points = 2500
    freqs = arange(num_points) * 1e10  # Hz
    freqs *= 1e-12  # THz
    
    times = arange(2 * (num_points - 1))
    times = times / (mean(diff(freqs)) * 2 * (num_points - 1))  # ps
    
    E_sim_w_abs = stats.lognorm.pdf(freqs, 0.5, scale=0.4)
    
    E_sim_w_abs /= max(E_sim_w_abs)
    
    # E_sim_w_abs = fromDb(30 + toDb(E_sim_w_abs))
    E_sim_w_arg = phase_factor(n_air, 0, 1.6e-3, freqs*1e12)
    # plot(freqs, E_sim_w_arg)
    E_sim_ref = irfft(E_sim_w_abs * E_sim_w_arg)  # * 100
    
    # for ns_floor in [-120, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10]:
    for ns_floor in [-90, -60, -30]:
        E_sim_ref += fromDb(ns_floor) * random.normal(0, 0.02, E_sim_ref.size)
        write_data(times, 100 * E_sim_ref, str(ns_floor) + '_ref', out_dir)  # THz


sim_refs()
