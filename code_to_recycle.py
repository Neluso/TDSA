
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

########################################################################################################################

# def fp_first_der_m(n, k, L, f, m):
#     n_cplx = n - 1j * k
#     rho_prima_rho_term = (2 * n_air) / (n_cplx ** 2 - n_air ** 2)
#     fp_val = 2 * m * (rho_prima_rho_term - 1j * 2 * pi * L / c_0)
#     fp_val *= fp_m(n, k, L, f, m)
#     # fp_val *= fp_full(n, k, L, f)
#     return fp_val
#
#
# def fp_second_der_m(f, m, n, k, L):
#     n_cplx = n - 1j * k
#     rho_prima_rho_term = (2 * n_air) / (n_cplx ** 2 - n_air ** 2)
#     fp_val = (2 * m * (rho_prima_rho_term - 1j * 2 * pi * f * L / c_0))**2
#     fp_val -= 8 * m * (4 * n_cplx * n_air / (n_cplx ** 2 - n_air ** 2))
#     fp_val *= fp_m(f, m, n, k, L)
#     return fp_val
########################################################################################################################

layers = 2
for i in range(1, layers - 1):
    print(layers - i - 1)
