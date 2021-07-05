from TDSA import *


def t_i_j(n_i, n_j):
    return (2 * n_i) / (n_i + n_j)


def r_i_j(n_i, n_j):
    return (n_i - n_j) / (n_i + n_j)


def phase_factor(n_i, thick_i, freq):  # theta in radians
    omg = 2 * pi * freq
    phi = omg * thick_i / c_0
    return exp(1j * n_i * phi)


def fabry_perot(n_i, thick_i, n_im1, n_ip1, freq):
    fp = 1 - r_i_j(n_im1, n_i) * r_i_j(n_i, n_ip1) * phase_factor(n_i, 2 * thick_i, freq)
    return 1 / fp


def H_sim(freq, n_i, n_im1, n_ip1, thick_i, H_prev):  # in first input: t_i_j(n_im1, n_i)
    H_i = H_prev * t_i_j(n_i, n_ip1) * phase_factor(n_i, thick_i, freq) * fabry_perot(n_i, thick_i, n_im1, n_ip1, freq)
    return H_i


def full_adiab_H_sim():
    return 0


file_path = './20210705_mesures/20210705_adiabatic_acrilic/'
sample_num = 1
t_ref, E_ref = read_1file(file_path + 'ref' + str(sample_num) + '.txt')
t_sam_ad, E_sam_ad = read_1file(file_path + 'adiab' + str(sample_num) + '.txt')
t_sam_na, E_sam_na = read_1file(file_path + 'noad' + str(sample_num) + '.txt')
delta_t_ref = mean(diff(t_ref))
enlargement = 10 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_ad = zero_padding(E_sam_ad, 0, enlargement)
t_sam_ad = concatenate((t_sam_ad, t_sam_ad[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_na = zero_padding(E_sam_na, 0, enlargement)
t_sam_na = concatenate((t_sam_na, t_sam_na[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
t_sam_ad *= 1e-12
t_sam_na *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam_ad, E_sam_ad_w = fourier_analysis(t_sam_ad, E_sam_ad)
f_sam_na, E_sam_na_w = fourier_analysis(t_sam_na, E_sam_na)


H_ad_w = E_sam_ad_w / E_ref_w
H_na_w = E_sam_na_w / E_ref_w



figure(1)
plot(t_ref * 1e12, E_ref, label='ref')
plot(t_sam_ad * 1e12, E_sam_ad, label='adiab')
plot(t_sam_na * 1e12, E_sam_na, label='no-ad')
xlim([0, 50])
legend()
savefig(file_path + str(sample_num) + 'traces.png')
figure(2)
plot(f_ref * 1e-12, abs(H_ad_w), label='adiab')
plot(f_ref * 1e-12, abs(H_na_w), label='no-ad')
xlim([0.05, 0.9])
ylim([0, 2])
legend()
savefig(file_path + str(sample_num) + 'transfer_func.png')
show()
