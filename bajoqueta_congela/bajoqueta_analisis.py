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


# t_ref, E_ref = read_1file('./data_bank/ref_seca.txt')
# t_sam, E_sam = read_1file('./data_bank/sam_seca.txt')
#
# t_ref *= 1e-12
# t_sam *= 1e-12
# f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
# n_baj_mean = zeros(f_sam.size)
# thickness = 0.47e-3  # m --- 470 um
# n_baj_aux, alpha_baj_aux, n_baj_avg_aux = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
#
# plot(f_sam * 1e-12, n_baj_mean)


for i in range(3):
    t_ref, E_ref = read_1file('./data_bank/ref_sec_2_' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/sam_sec_2_' + str(i + 1) + '.txt')
    
    t_ref *= 1e-12
    t_sam *= 1e-12
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        n_baj_mean = zeros(f_sam.size)
        alpha_baj_mean = zeros(f_sam.size)
    thickness = 0.47e-3  # m --- 470 um
    n_baj_aux, alpha_baj_aux, n_baj_avg_aux = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
    n_baj_mean += n_baj_aux
    alpha_baj_mean += alpha_baj_aux
figure(1)
plot(f_sam * 1e-12, n_baj_mean / 3, label='24h', lw=0.5)
figure(2)
plot(f_sam * 1e-12, 1e-2 * alpha_baj_mean / 3, label='24h', lw=0.5)


for i in range(3):
    t_ref, E_ref = read_1file('./data_bank/ref_sec_3_' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/sam_sec_3_' + str(i + 1) + '.txt')
    t_ref *= 1e-12
    t_sam *= 1e-12
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        n_baj_mean = zeros(f_sam.size)
        alpha_baj_mean = zeros(f_sam.size)
    thickness = 0.36e-3  # m --- 360 um
    n_baj_aux, alpha_baj_aux, n_baj_avg_aux = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
    n_baj_mean += n_baj_aux
    alpha_baj_mean += alpha_baj_aux
figure(1)
plot(f_sam * 1e-12, n_baj_mean / 3, label='48h', lw=0.5)
figure(2)
plot(f_sam * 1e-12, 1e-2 * alpha_baj_mean / 3, label='48h', lw=0.5)

figure(1)
xlim([0.1, 1.6])
ylim([0.9, 2.1])
xlabel(r'$f\ (THz)$')
ylabel(r'n')
legend()
savefig('./n_baj.png')
clf()
figure(2)
xlim([0.1, 1.6])
ylim([0, 200])
xlabel(r'$f\ (THz)$')
ylabel(r'$\alpha\ (cm^{-1})$')
legend()
savefig('./alpha_baj.png')
clf()

for i in range(7):
    t_ref, E_ref = read_1file('./data_bank/1a_descongelació/ref' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/1a_descongelació/sam' + str(i + 1) + '.txt')
    t_ref *= 1e-12
    t_sam *= 1e-12
    # E_sam = zero_padding(E_sam, 0, E_sam.size)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        E_sam_w_max = max(abs(E_sam_w))
    n_1a, alpha_1a, n_avg_1a = jepsen_index(t_ref, E_ref, t_sam, E_sam, 3.3e-3)
    figure(1)
    plot(f_sam * 1e-12, n_1a, label=str(i + 1), lw=0.5)
    figure(2)
    plot(f_sam * 1e-12, alpha_1a * 1e-2, label=str(i + 1), lw=0.5)
figure(1)
xlim([0.1, 1.2])
ylim([0.9, 3.5])
xlabel(r'$f\ (THz)$')
ylabel(r'n')
legend()
savefig('./n_unfreeze_1.png')
clf()
figure(2)
xlim([0.1, 1.2])
ylim([0, 80])
xlabel(r'$f\ (THz)$')
ylabel(r'$\alpha\ (cm^{-1})$')
legend()
savefig('./alpha_unfreeze_1.png')
clf()


for i in range(7, 10):
    t_ref, E_ref = read_1file('./data_bank/2a_descongelació/ref' + str(i + 1) + '.txt')
    t_sam, E_sam = read_1file('./data_bank/2a_descongelació/sam' + str(i + 1) + '.txt')
    t_ref *= 1e-12
    t_sam *= 1e-12
    # E_sam = zero_padding(E_sam, 0, E_sam.size)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    if i == 0:
        E_sam_w_max = max(abs(E_sam_w))
    n_1a, alpha_1a, n_avg_1a = jepsen_index(t_ref, E_ref, t_sam, E_sam, 3.3e-3)
    figure(3)
    plot(f_sam * 1e-12, n_1a, label=str(i + 1), lw=0.5)
    figure(4)
    plot(f_sam * 1e-12, alpha_1a * 1e-2, label=str(i + 1), lw=0.5)
figure(3)
xlim([0.1, 1.2])
ylim([0.9, 3.5])
xlabel(r'$f\ (THz)$')
ylabel(r'n')
legend()
savefig('./n_unfreeze_2.png')
clf()
figure(4)
xlim([0.1, 1.2])
ylim([0, 80])
xlabel(r'$f\ (THz)$')
ylabel(r'$\alpha\ (cm^{-1})$')
legend()
savefig('./alpha_unfreeze_2.png')
clf()


time_points = list()
E_sam_w_power = list()
for i in range(30):
    try:
        t_sam, E_sam = read_1file('./data_bank/descong_invivo/' + str(i + 1) + '.txt')
    except:
        continue
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    time_points.append(i)
    E_sam_w_power.append(abs(dot(E_sam_w, conjugate(E_sam_w))))
time_points = array(time_points)
E_sam_w_power = array(E_sam_w_power)
figure(1)
plot(time_points, toDb_0(E_sam_w_power), lw=0.5)
xlabel(r'$\Delta t\ (min)$')
ylabel(r'$T\ (dB)$')
savefig('./evol_T.png')
clf()
