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
sample_num = 0
for sample_num in range(20):
    t_ref, E_ref = read_1file(file_path + 'ref.txt')
    t_sam, E_sam = read_1file(file_path + str(sample_num) + '.txt')
    delta_t_ref = mean(diff(t_ref))
    enlargement = 10 * E_ref.size
    E_ref = zero_padding(E_ref, 0, enlargement)
    t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    E_sam = zero_padding(E_sam, 0, enlargement)
    t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    t_ref *= 1e-12
    t_sam *= 1e-12


    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)


    H_w = E_sam_w / E_ref_w

    # H_teo = H_sim(f_ref, 1.45, n_air, 1.50, 500e-6, t_i_j(n_air, 1.45))
    # H_teo = H_sim(f_ref, 1.475, 1.45, 1.50, 500e-6, H_teo)
    # H_teo = H_sim(f_ref, 1.55, 1.50, n_air, 500e-6, H_teo)


    # figure(1)
    # plot(t_ref * 1e12, E_ref, label='ref')
    # plot(t_sam * 1e12, E_sam, label='adiab')
    # plot(t_ref * 1e12, irfft(H_teo * E_ref_w), label='sim')
    # xlim([0, 50])
    # legend()
    # savefig(file_path + str(sample_num) + 'traces.png')
    figure(2)
    plot(f_ref * 1e-12, abs(H_w), label='adiab')
    # plot(f_ref * 1e-12, abs(H_teo), label='sim')
    xlim([0.05, 0.9])
    ylim([0, 2])
    # legend()
    # # savefig(file_path + str(sample_num) + 'transfer_func.png')
    # figure(3)
    # plot(sample_num, sum(abs(E_ref) - abs(E_sam)), 'b.')
    # ylim([0, 60])
show()
