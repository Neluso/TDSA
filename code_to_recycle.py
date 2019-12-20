
# k = res.x
# pulses = list()
# for i in range(int(round(k.size / 2))):
#     pulses.append((k[2 * i + 1], k[2 * i]))
# pulses = sorted(pulses)
# delta_t = abs(pulses[0][0] - pulses[1][0])
# thickness = c_0 * delta_t * 1e-12 / (2 * 1.51 * cos(pi / 6))  # m
# thickness *= 1e3  # mm
# print(thickness, 'mm')
#
# plot(t_sam, E_sam, lw=1)
# E_sim = objective_function(k, E_ref, E_sam, t_ref)
# plot(t_ref, E_sim, lw=1)
# show()


########################################################################################################################


# def objective_function(k, *args):
#     k0, k1, k2, k3, k4, k5, k6, k7 = k
#     E_val, E_meas, t_val = args
#     idx_k1 = centre_loc(E_val) - int(round(k1 / mean(diff(t_val))))
#     idx_k3 = centre_loc(E_val) - int(round(k3 / mean(diff(t_val))))
#     idx_k5 = centre_loc(E_val) - int(round(k5 / mean(diff(t_val))))
#     idx_k7 = centre_loc(E_val) - int(round(k7 / mean(diff(t_val))))
#     return k0 * roll(E_val, - idx_k1) + k2 * roll(E_val, - idx_k3) + k4 * roll(E_val, - idx_k5) + k6 * roll(E_val, - idx_k7)


layers = 2
for i in range(1, layers - 1):
    print(layers - i - 1)
