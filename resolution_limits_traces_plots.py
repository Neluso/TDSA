from TDSA import *
import os


in_dir = './output/traces/'
dir_list = os.listdir(in_dir)

for file_i in dir_list:
    t_sim, E_sim = read_1file(in_dir + file_i)  # t_ref in s
    f_sim, E_sim_w = fourier_analysis(t_sim, E_sim)  # f_ref in Hz
    plot(t_sim, E_sim)

show()
