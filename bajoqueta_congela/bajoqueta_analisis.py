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


for i in range(3):
    t_ref, E_ref = read_1file('./data_bank/ref_sec_2_' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/sam_sec_2_' + str(i + 1) + '.txt')
    
    t_ref *= 1e-12
    t_sam *= 1e-12
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        n_baj_mean = zeros(f_sam.size)
    thickness = 0.47e-3  # m --- 470 um
    n_baj_aux, alpha_baj_aux, n_baj_avg_aux = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
    n_baj_mean += n_baj_aux
plot(f_sam * 1e-12, n_baj_mean / 3, 'b')


for i in range(3):
    t_ref, E_ref = read_1file('./data_bank/ref_sec_3_' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/sam_sec_3_' + str(i + 1) + '.txt')
    t_ref *= 1e-12
    t_sam *= 1e-12
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        n_baj_mean = zeros(f_sam.size)
    thickness = 0.36e-3  # m --- 360 um
    n_baj_aux, alpha_baj_aux, n_baj_avg_aux = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
    n_baj_mean += n_baj_aux
plot(f_sam * 1e-12, n_baj_mean / 3, 'r')
xlim([0.1, 1.6])
ylim([0.9, 2.1])
# for i in range(7):
    # t_ref, E_ref = read_1file('./data_bank/1a_descongelació/ref' + str(i + 1) + '.txt')
    # t_sam, E_sam = read_1file('./data_bank/1a_descongelació/sam' + str(i + 1) + '.txt')
    # E_sam = zero_padding(E_sam, 0, E_sam.size)
    # f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    # if i == 0:
    #     E_sam_w_max = max(abs(E_sam_w))
#     figure(2)
#     plot(f_sam, prettyfy(E_sam_w, E_sam_w_max), lw=0.5, label=str(i + 1))
    # plot(f_sam, self_unwrap(t_sam, E_sam))  #, label=str(i + 1))

# time_points = list()
# E_sam_w_power = list()
# for i in range(30):
#     try:
#         t_sam, E_sam = read_1file('./data/' + str(i + 1) + '.txt')
#     except:
#         continue
#     E_sam = zero_padding(E_sam, 0, E_sam.size)
#     f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
#     if i == 0:
#         E_sam_w_max = max(abs(E_sam_w))
#     # if i == 0:
#     #     E_sam_w_max = max(abs(E_sam_w))
#     # figure(1)
#     # print(dot(E_sam_w, conjugate(E_sam_w)))
#     # quit()
#     time_points.append(i)
#     E_sam_w_power.append(abs(dot(E_sam_w / E_sam_w_max, conjugate(E_sam_w))))
# time_points = array(time_points)
# E_sam_w_power = array(E_sam_w_power)
# plot(time_points, toDb(E_sam_w_power), lw=0.5)
# xlabel('mins')
# ylabel('P (dB)')
#
#
# # legend()
show()
