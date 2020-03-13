from TDSA import *
from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
from matplotlib.pyplot import *
from scipy.signal import windows
from scipy.optimize import differential_evolution
from scipy.special import erf
from scipy.interpolate import Akima1DInterpolator


particle_count = 150
max_iterations = 1500


wh = open('./output/stability_results.csv', 'w')
# wh_test = open('./output/notch_results_test.csv', 'w')
all_r = list()


for i in range(10):
    t_ref, E_ref = read_slow_data('./data/notch/camp1/ref' + str(i + 1) + '.txt')
    # t_ref_test, E_ref_test = read_slow_data('./data/notch/ref_20200309.txt')
    t_ref_test, E_ref_test = read_slow_data('./data/notch/ref_20200310_1.txt')
    win_idx = int(round(E_ref.size / 3))
    window = windows.tukey(win_idx, 0.2)
    window = zero_padding(window, 0, E_ref.size - win_idx)
    
    E_ref *= window
    # E_pap *= window
    E_ref_test *= window
    
    delta_t_ref = mean(diff(t_ref))
    enlargement = 2 * E_ref.size

    # plot(t_ref, E_ref)
    # plot(t_ref_test, E_ref_test)
    # plot(arange(window.size) * delta_t_ref, window)
    # show()
    # quit()
    
    E_ref = zero_padding(E_ref, 0, enlargement)
    t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    E_ref_test = zero_padding(E_ref_test, 0, enlargement)
    t_ref_test = concatenate((t_ref_test, t_ref_test[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))

    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_ref_test, E_ref_test_w = fourier_analysis(t_ref_test, E_ref_test)

    # f_min, f_max = f_min_max_idx(f_ref * 1e12, 0.2, 1.0)
    # f_ref = f_ref[f_min:f_max]
    # E_ref_w = E_ref_w[f_min:f_max]
    # E_ref_test_w = E_ref_test_w[f_min:f_max]
    
    all_r.append(E_ref_w/E_ref_test_w)
    

all_r = array(all_r)

errorbar(f_ref, mean(all_r, 0), yerr=std(all_r, 0))
show()
