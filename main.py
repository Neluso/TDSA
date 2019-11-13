from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, title
from read_data import read_data
from jepsen_index import jepsen_index
from numpy import array


t_ref, E_ref = read_data('ref.txt')
t_sam, E_sam = read_data('pp 0.5cnt 1.95mm.txt')
t_ref *= 1e-12  # seconds
t_sam *= 1e-12  # seconds
thickness = 1.95e-3  # m
f_ref, f_min_idx, f_max_idx, n, delta_phi = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
bins = int(f_ref.size/8)


print('Plotting data')
figure(1)
plot(f_ref[f_min_idx:2*f_max_idx], n[f_min_idx:2*f_max_idx])
title('Refractive index')
figure(2)
plot(f_ref[f_min_idx:2*f_max_idx], delta_phi[f_min_idx:2*f_max_idx])
title('delta_phi')
show()
print('Closing')