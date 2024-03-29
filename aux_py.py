from numpy import *
from matplotlib.pyplot import *


ref_pp_t = 7.59
theta_t = array((0, 10, 20, 22.5, 30, 40, 45, 50, 60, 67.5, 70))  # º
theta_r = array((30, 35, 40, 45, 50, 55, 60, 65))  # º
v_pp_t = array((5.52, 5.53, 5.43, 5.38, 5.23, 4.96, 4.81, 4.59, 2.72, 1.11, 1.06))  # V
ref_pp_r = array((4.24, 4.49, 4.23, 4.37, 4.33, 4.36, 4.29, 4.11))
v_pp_r = array((1.94, 1.93, 1.66, 1.53, 1.30, 1.07, 0.87, 0.88))
v_pp_err = 0.01  # * ones(v_pp.size)  # V
theta_err = 2  # * ones(theta.size)  # º

theta_total = array((30, 40, 45, 50, 60))
v_pp_total = list()
v_pp_total_err = list()

for i in range(theta_t.size):
    for j in range(theta_r.size):
        if theta_t[i] == theta_r[j]:
            v_pp_total.append(v_pp_t[i] / ref_pp_t + v_pp_r[j] / ref_pp_r[j])
            v_pp_total_err.append(sqrt((v_pp_err / ref_pp_t)**2 + (v_pp_err / ref_pp_r[j])**2))
v_pp_total = array(v_pp_total)
v_pp_total_err = array(v_pp_total_err)


figure(1)
errorbar(theta_t, v_pp_t / ref_pp_t, yerr=v_pp_err / ref_pp_t, xerr=0*theta_err, lw=1, label='t')
errorbar(theta_r, v_pp_r / ref_pp_r, yerr=v_pp_err / ref_pp_r, xerr=0*theta_err, lw=1, label='r')
errorbar(theta_total, v_pp_total, yerr=v_pp_total_err, xerr=0*theta_err, lw=1, label='sum')
axvline(55, color='r', linestyle='dashed', lw=0.5)
hlines(1, 0, 90, colors='r', linestyles='dashed', lw=0.5)
ylim([0, 1.5])
xlim([0, 70])
legend()
show()

# ref_pulse_idx = centre_loc(E_ref)  TODO review this part of code
# # lower bounds for k0, k1, k2, k3 respectively
# k_min = [0, t_ref[0], -1, t_ref[ref_pulse_idx]]
# # upper bounds for k0, k1, k2, k3 respectively
# k_max = [1, t_ref[ref_pulse_idx], 1, t_ref[-1]]
#
# k0 = list()
# k1 = list()
# k2 = list()
# k3 = list()
#
# thick = list()
# error_func = list()
# repetitions = 10
#
# for i in range(repetitions):
#
#     print('Iteration', i + 1, 'of', repetitions)
#
#     k, fopt = pso(min_function, k_min, k_max,
#                   args=(E_ref, E_sam, t_ref),
#                   swarmsize=1000,
#                   maxiter=2000,
#                   f_ieqcons=constraints,
#                   phig=0.1,
#                   phip=0.1,
#                   minstep=1e-10,
#                   minfunc=1e-10,
#                   debug=False)
#
#
#     # figure(1)
#     # title('Fit')
#     # plot(t_sam, E_sam, lw=1, label='sam')
#     # plot(t_ref, objective_function(k, *(E_ref, E_sam, t_ref)), lw=1, label='fit')
#     # legend()
#
#     delta_t = k[3] - k[1]
#     thickness = c_0 * delta_t * 1e-12 / (2 * 2.6)  # m
#     thickness *= 1e3  # mm
#     thick.append(thickness)
#     error_func.append(fopt)
#     k0.append(k[0])
#     k1.append(k[1])
#     k2.append(k[2])
#     k3.append(k[3])
#     # print('k0 =', k[0])
#     # print('k1 =', k[1])
#     # print('k2 =', k[2])
#     # print('k3 =', k[3])
#     # print('delta t =', delta_t, 'ps')
#     # print('d =', thickness, 'mm')
#
# figure(1)
# plot(arange(repetitions), array(thick))
# figure(2)
# plot(arange(repetitions), array(error_func))
# show()
