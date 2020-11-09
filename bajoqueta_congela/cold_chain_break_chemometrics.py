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


for i in range(5):
    thick = open('./cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/thickness.txt')
    thick = float(thick.read())
    thick = 1e3  # mm
    for j in range(4):
        try:
            t_ref, E_ref = read_1file(
                './cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/rotura_' + str(j + 1) + '/ref1.txt')
            t_sam, E_sam = read_1file(
                './cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/rotura_' + str(j + 1) + '/sam1.txt')
        except:
            continue
        
        t_ref *= 1e-12
        t_sam *= 1e-12

        f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

        H_w = E_sam_w / E_ref_w

        plot(f_ref, gradient(gradient(abs(H_w))))

        show()
