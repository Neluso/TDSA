from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, title
from read_data import read_data
from jepsen_index import jepsen_index


t_ref, E_ref = read_data('ref.txt')
t_sam, E_sam = read_data('pp 0.5cnt 1.95mm.txt')
t_ref *= 1e-12  # seconds
t_sam *= 1e-12  # seconds
thickness = 1.95e-3  # m
f_ref, f_min_idx, f_max_idx, n, alpha_f = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
bins = int(f_ref.size/8)


print('Plotting data')
figure(1)
plot(f_ref[f_min_idx:2*f_max_idx], 0.01 * alpha_f[f_min_idx:2*f_max_idx])
title('Absorption')
xlabel('f (THz)')
ylabel(r'$\alpha (cm^{-1})$')
legend()
show()
print('Closing')