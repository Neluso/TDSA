from TDSA import *


t_ref_1, E_ref_1 = read_1file('./data/er_caba/v2/ref1_stitch_new.txt')
t_sam_1, E_sam_1 = read_1file('./data/er_caba/v2/sam1_stitch_new.txt')
t_ref_2, E_ref_2 = read_1file('./data/er_caba/v2/ref2_stitch_new.txt')
t_sam_2, E_sam_2 = read_1file('./data/er_caba/v2/sam2_stitch_new.txt')
t_ref_3, E_ref_3 = read_1file('./data/er_caba/v2/ref3_stitch_new.txt')
t_sam_3, E_sam_3 = read_1file('./data/er_caba/v2/sam3_stitch_new.txt')


delta_t_ref = mean(diff(t_ref_1))
enlargement = 10 * E_ref_1.size
E_ref_1 = zero_padding(E_ref_1, 0, enlargement)
t_ref_1 = concatenate((t_ref_1, t_ref_1[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_1 = zero_padding(E_sam_1, 0, enlargement)
t_sam_1 = concatenate((t_sam_1, t_sam_1[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))

delta_t_ref = mean(diff(t_ref_2))
enlargement = 10 * E_ref_3.size
E_ref_2 = zero_padding(E_ref_2, 0, enlargement)
t_ref_2 = concatenate((t_ref_2, t_ref_2[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_2 = zero_padding(E_sam_2, 0, enlargement)
t_sam_2 = concatenate((t_sam_2, t_sam_2[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))

delta_t_ref = mean(diff(t_ref_3))
enlargement = 10 * E_ref_3.size
E_ref_3 = zero_padding(E_ref_3, 0, enlargement)
t_ref_3 = concatenate((t_ref_3, t_ref_1[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_3 = zero_padding(E_sam_3, 0, enlargement)
t_sam_3 = concatenate((t_sam_3, t_sam_3[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))


# plot(t_ref, E_ref)
# plot(t_sam, E_sam)
# show()
# quit()
t_ref_1 *= 1e-12
t_sam_1 *= 1e-12
f_ref_1, E_ref_w_1 = fourier_analysis(t_ref_1, E_ref_1)
f_sam_1, E_sam_w_1 = fourier_analysis(t_sam_1, E_sam_1)

t_ref_2 *= 1e-12
t_sam_2 *= 1e-12
f_ref_2, E_ref_w_2 = fourier_analysis(t_ref_2, E_ref_2)
f_sam_2, E_sam_w_2 = fourier_analysis(t_sam_2, E_sam_2)

t_ref_3 *= 1e-12
t_sam_3 *= 1e-12
f_ref_3, E_ref_w_3 = fourier_analysis(t_ref_3, E_ref_3)
f_sam_3, E_sam_w_3 = fourier_analysis(t_sam_3, E_sam_3)

figure(1)
plot(f_ref_1, toDb(E_ref_w_1), 'b')
plot(f_sam_1, toDb(E_sam_w_1), 'g')
figure(2)
plot(f_ref_2, toDb(E_ref_w_2), 'b')
plot(f_sam_2, toDb(E_sam_w_2), 'g')
figure(3)
plot(f_ref_3, toDb(E_ref_w_3), 'b')
plot(f_sam_3, toDb(E_sam_w_3), 'g')
show()
quit()

f_min_idx, f_max_idx = f_min_max_idx(f_ref_1, 0, 5)

absH_w_1 = abs(E_sam_w_1[f_min_idx:f_max_idx] / E_ref_w_1[f_min_idx:f_max_idx])
angH_w_1 = jepsen_unwrap(t_ref_1, E_ref_1, t_sam_1, E_sam_1)

absH_w_2 = abs(E_sam_w_2[f_min_idx:f_max_idx] / E_ref_w_2[f_min_idx:f_max_idx])
angH_w_2 = jepsen_unwrap(t_ref_2, E_ref_2, t_sam_2, E_sam_2)

plot(f_ref_1[f_min_idx:f_max_idx], absH_w_1)
plot(f_ref_2[f_min_idx:f_max_idx], absH_w_2)
# plot(f_ref_3[f_min_idx:f_max_idx], abs(E_ref_w_2 / E_ref_w_2))
show()

# f_min_idx, f_max_idx = f_min_max_idx(f_ref_1, 0.1, 0.4)
# f_ref_1 *= 1e-12
# f_sam_1 *= 1e-12
#
# figure(1)
# if i == 0:
#     plot(f_ref, prettyfy(E_ref_w, amax(abs(E_ref_w))), label='ref')
# plot(f_sam, prettyfy(E_sam_w, amax(abs(E_ref_w))), label='sam')
# xlabel('f (THz)')
# legend()
#
# figure(2)
# plot(f_ref[f_min_idx:f_max_idx], absH_w[f_min_idx:f_max_idx])
# ylim([0, 1.5])
# xlabel('f (THz)')
#
# figure(3)
# plot(f_ref[f_min_idx:f_max_idx], angH_w[f_min_idx:f_max_idx])
# xlabel('f (THz)')
#
# show()
