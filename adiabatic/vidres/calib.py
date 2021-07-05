from TDSA import *
from jepsen_index import *


sample_num = 4
t_ref, E_ref = read_1file('./20210702_adiab_vidres/calib_adiab/ref' + str(sample_num) + '.txt')
t_sam_r, E_sam_r = read_1file('./20210702_adiab_vidres/calib_adiab/roig' + str(sample_num) + '.txt')
t_sam_v, E_sam_v = read_1file('./20210702_adiab_vidres/calib_adiab/verd' + str(sample_num) + '.txt')
delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_r = zero_padding(E_sam_r, 0, enlargement)
t_sam_r = concatenate((t_sam_r, t_sam_r[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_v = zero_padding(E_sam_v, 0, enlargement)
t_sam_v = concatenate((t_sam_v, t_sam_v[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
t_sam_r *= 1e-12
t_sam_v *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam_r, E_sam_r_w = fourier_analysis(t_sam_r, E_sam_r)
f_sam_v, E_sam_v_w = fourier_analysis(t_sam_v, E_sam_v)


H_r_w = E_sam_r_w / E_ref_w
H_v_w = E_sam_v_w / E_ref_w

print(n_avg_calc(t_ref, E_ref, t_sam_r, E_sam_r, 0.04e-3))
print(n_avg_calc(t_ref, E_ref, t_sam_v, E_sam_v, 0.03e-3))
plot(t_ref, E_ref, label='ref')
plot(t_sam_r, E_sam_r, label='roig')
plot(t_sam_v, E_sam_v, label='verd')
legend()


# plot(jepsen_unwrap(t_ref, E_ref, t_sam_r, E_sam_r))
# plot(jepsen_unwrap(t_ref, E_ref, t_sam_v, E_sam_v))
show()

