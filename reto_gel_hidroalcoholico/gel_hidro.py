from TDSA import *


t_ref, E_ref = read_1file('./ref.txt')
t_sam, E_sam = read_1file('./sam.txt')

t_ref *= 1e-12
t_sam *= 1e-12

n_gel, alpha_gel, n_gel_avg = jepsen_index(t_sam, E_ref, t_sam, E_sam, 90*1e-6, 2.6)

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
figure(1)
plot(f_ref * 1e-12, n_gel)
figure(2)
plot(f_ref * 1e-12, alpha_gel)
show()
