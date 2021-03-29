from TDSA import *


t_ref, E_ref = read_1file('ref1.txt')
t_sam, E_sam = read_1file('PP1.txt')
t_ref *= 1e-12
t_sam *= 1e-12
enlargement = 9 * E_ref.size
delta_t_ref = mean(diff(t_ref))
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
n, alpha, n_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, 1.92e-3)  # 1.94 mm

# plot(f_ref, abs(alpha)*1e-2)

f_PP, alpha_PP = read_1file('alpha_PP.csv')
f_PP *= 1e12

plot(f_PP, 0.6*(alpha_PP-1))
wh = open('alpha_PP_mod.csv', 'w')
for i in range(f_PP.size):
    wh.write(str(f_PP[i]*1e-12) + ',' + str(0.6*(alpha_PP[i]-1)) + '\n')

show()
