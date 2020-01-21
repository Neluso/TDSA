from read_data import *
from DSP_functions import *
from aux_functions import *


t_ref_10_01, E_ref_10_01 = read_slow_data('./data/old/cork/ref10_01.txt')


t_ref_10_05, E_ref_10_05 = read_slow_data('./data/old/cork/ref10_05.txt')
t_ref30_05, E_ref30_05 = read_slow_data('./data/old/cork/ref30_05.txt')
t_ref_100_05, E_ref_100_05 = read_slow_data('./data/old/cork/ref100_05.txt')


f_ref_10_01, E_ref_w10_01 = fourier_analysis(t_ref_10_01, E_ref_10_01, E_ref_10_01.size)

f_ref_10_05, E_ref_w10_05 = fourier_analysis(t_ref_10_05, E_ref_10_05, E_ref_10_05.size)
f_ref30_05, E_ref_w30_05 = fourier_analysis(t_ref30_05, E_ref30_05, E_ref30_05.size)
f_ref_100_05, E_ref_w100_05 = fourier_analysis(t_ref_100_05, E_ref_100_05, E_ref_100_05.size)


figure(1)

plot(f_ref_10_05, toDb(E_ref_w10_05), lw=0.8, label='10_05')
plot(f_ref30_05, toDb(E_ref_w30_05), lw=0.8, label='30_05')
plot(f_ref_100_05, toDb(E_ref_w100_05), lw=0.8, label='100_05')
plot(f_ref_10_01, toDb(E_ref_w10_01), lw=0.8, label='10_01')
legend()

show()
