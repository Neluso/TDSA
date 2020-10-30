from TDSA import *


def centroid_E2(t_val, E_val):  # t_centroid
    t_centroid = sum(t_val * abs(E_val)**2) / sum(abs(E_val)**2)
    t_idx = where(t_val <= t_centroid)[0]
    return t_idx[-1]


def jepsen_unwrap(t_ref, E_ref, t_sam, E_sam):
    t_ref_0 = centroid_E2(t_ref, E_ref)
    t_sam_0 = centroid_E2(t_sam, E_sam)
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    phi_ref_0 = 2 * pi * f_ref * t_ref_0
    phi_sam_0 = 2 * pi * f_sam * t_sam_0
    phi_ref_0_red = angle(E_ref_w * exp(-1j * phi_ref_0))
    phi_sam_0_red = angle(E_sam_w * exp(-1j * phi_sam_0))
    return unwrap(phi_sam_0_red - phi_ref_0_red)


def self_unwrap(t_sam, E_sam):
    t_sam_0 = centroid_E2(t_sam, E_sam)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    phi_sam_0 = 2 * pi * f_sam * t_sam_0
    phi_sam_0_red = angle(E_sam_w * exp(-1j * phi_sam_0))
    return unwrap(phi_sam_0_red)


for i in range(4):
    thick = open('./cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/thickness.txt')
    thick = float(thick.read())
    # print(type(thick))
    # quit()
    for j in range(2):
        t_ref, E_ref = read_1file(
            './cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/rotura_' + str(j + 1) + '/ref1.txt')
        t_sam, E_sam = read_1file(
            './cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/rotura_' + str(j + 1) + '/sam1.txt')
        t_ref *= 1e-12
        t_sam *= 1e-12

        # delta_t_ref = mean(diff(t_ref))
        # enlargement = 9 * E_ref.size
        # E_ref = zero_padding(E_ref, 0, enlargement)
        # t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
        # E_sam = zero_padding(E_sam, 0, enlargement)
        # t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
        
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
        n_1a, alpha_1a, n_avg_1a = jepsen_index(t_ref, E_ref, t_sam, E_sam, thick*1e-3)
        figure(1 * (i+1))
        plot(f_sam * 1e-12, n_1a, label=str(j + 1) + r' $n_{avg}=$' + str(n_avg_1a), lw=0.5)
        xlim([0.1, 1.5])
        ylim([1, 3])
        legend()
        # figure(5 * (i+1))
        # plot(f_sam * 1e-12, alpha_1a * 1e-2, label=str(i + 1), lw=0.5)
show()
