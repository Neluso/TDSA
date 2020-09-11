# import numpy as np
# import read_data as rd
# import matplotlib.pyplot as plt
# import jepsen_index as jpind
# import DSP_functions as DSPf
import TDSA


muestra_j = '1'

for i in range(3):
    t_ref, E_ref = read_1file('./data/muestra_' + str(muestra_j) + '/ref_' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data/muestra_' + str(muestra_j) + '/sam_' + str(i + 1) + '.txt')
    # n_idx, alpha_r, n_idx_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, 0.75)
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    figure(1)
    if i == 0:
        plot(t_ref, E_ref)
    plot(t_sam, E_sam)
    figure(2)
    if i == 0:
        plot(f_ref, toDb(E_ref_w))
    plot(f_sam, toDb(E_sam_w))
show()
