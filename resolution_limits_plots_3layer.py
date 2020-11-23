from TDSA import *
import csv


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


# fh = open('resolution_limit.csv', 'r')
# csvfile = csv.reader('resolution_limit.csv')

d_mat_i = list()
d_mat_mean_i = list()
d_mat_std_i = list()

e_s_i = list()
e_s_mean_i = list()
e_s_std_i = list()

e_inf_i = list()
e_inf_mean_i = list()
e_inf_std_i = list()

tau_i = list()
tau_mean_i = list()
tau_std_i = list()

d_mat_m = list()
d_mat_mean_m = list()
d_mat_std_m = list()

e_s_m = list()
e_s_mean_m = list()
e_s_std_m = list()

e_inf_m = list()
e_inf_mean_m = list()
e_inf_std_m = list()

tau_m = list()
tau_mean_m = list()
tau_std_m = list()

d_mat_o = list()
d_mat_mean_o = list()
d_mat_std_o = list()

e_s_o = list()
e_s_mean_o = list()
e_s_std_o = list()

e_inf_o = list()
e_inf_mean_o = list()
e_inf_std_o = list()

tau_o = list()
tau_mean_o = list()
tau_std_o = list()

d_air = list()
d_air_mean = list()
d_air_std = list()

f_cutoff = list()
pDr = list()

# rows = list()

sample_dir = 'test_0'
error_dir = '5%'

sample_dir = sample_dir + '/' + error_dir

with open('./output/simulation_results/' + sample_dir + '/resolution_limit.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        
        d_mat_i.append(float(row[0]))
        d_mat_mean_i.append(float(row[1]))
        d_mat_std_i.append(float(row[2]))
        e_s_i.append(float(row[3]))
        e_s_mean_i.append(float(row[4]))
        e_s_std_i.append(float(row[5]))
        e_inf_i.append(float(row[6]))
        e_inf_mean_i.append(float(row[7]))
        e_inf_std_i.append(float(row[8]))
        tau_i.append(float(row[9]))
        tau_mean_i.append(float(row[10]))
        tau_std_i.append(float(row[11]))
        d_mat_m.append(float(row[12]))
        d_mat_mean_m.append(float(row[13]))
        d_mat_std_m.append(float(row[14]))
        e_s_m.append(float(row[15]))
        e_s_mean_m.append(float(row[16]))
        e_s_std_m.append(float(row[17]))
        e_inf_m.append(float(row[18]))
        e_inf_mean_m.append(float(row[19]))
        e_inf_std_m.append(float(row[20]))
        tau_m.append(float(row[21]))
        tau_mean_m.append(float(row[22]))
        tau_std_m.append(float(row[23]))
        d_mat_o.append(float(row[24]))
        d_mat_mean_o.append(float(row[25]))
        d_mat_std_o.append(float(row[26]))
        e_s_o.append(float(row[27]))
        e_s_mean_o.append(float(row[28]))
        e_s_std_o.append(float(row[29]))
        e_inf_o.append(float(row[30]))
        e_inf_mean_o.append(float(row[31]))
        e_inf_std_o.append(float(row[32]))
        tau_o.append(float(row[33]))
        tau_mean_o.append(float(row[34]))
        tau_std_o.append(float(row[35]))
        d_air.append(float(row[36]))
        d_air_mean.append(float(row[37]))
        d_air_std.append(float(row[38]))
        f_cutoff.append(float(row[39]))
        pDr.append(float(row[40]))


# print(pDr)
# quit()
d_mat_i = array(d_mat_i)  # * 1e6
d_mat_mean_i = array(d_mat_mean_i)  # * 1e6
d_mat_std_i = array(d_mat_std_i)  # * 1e6
d_mat_pDr_i = array(pDr)
e_s_i = array(e_s_i)
e_s_mean_i = array(e_s_mean_i)
e_s_std_i = array(e_s_std_i)
e_inf_i = array(e_inf_i)
e_inf_mean_i = array(e_inf_mean_i)
e_inf_std_i = array(e_inf_std_i)
tau_i = array(tau_i)
tau_mean_i = array(tau_mean_i)
tau_std_i = array(tau_std_i)
d_mat_m = array(d_mat_m)  # * 1e6
d_mat_mean_m = array(d_mat_mean_m)  # * 1e6
d_mat_std_m = array(d_mat_std_m)  # * 1e6
d_mat_pDr_m = array(pDr)
e_s_m = array(e_s_m)
e_s_mean_m = array(e_s_mean_m)
e_s_std_m = array(e_s_std_m)
e_inf_m = array(e_inf_m)
e_inf_mean_m = array(e_inf_mean_m)
e_inf_std_m = array(e_inf_std_m)
tau_m = array(tau_m)
tau_mean_m = array(tau_mean_m)
tau_std_m = array(tau_std_m)
d_mat_o = array(d_mat_o)  # * 1e6
d_mat_mean_o = array(d_mat_mean_o)  # * 1e6
d_mat_std_o = array(d_mat_std_o)  # * 1e6
d_mat_pDr_o = array(pDr)
e_s_o = array(e_s_o)
e_s_mean_o = array(e_s_mean_o)
e_s_std_o = array(e_s_std_o)
e_inf_o = array(e_inf_o)
e_inf_mean_o = array(e_inf_mean_o)
e_inf_std_o = array(e_inf_std_o)
tau_o = array(tau_o)
tau_mean_o = array(tau_mean_o)
tau_std_o = array(tau_std_o)

d_air = array(d_air)
d_air_mean = array(d_air_mean)
d_air_std = array(d_air_std)
f_cutoff = array(f_cutoff)
pDr = array(pDr)
lmda = (c_0 / (f_cutoff * 1e12)) * 1e6  # um
lmda10 = list()
lmda20 = list()
lmda30 = list()
lmda40 = list()
lmda50 = list()
lmda60 = list()
lmda90 = list()

d_mat_mean10 = list()
d_mat_mean20 = list()
d_mat_mean30 = list()
d_mat_mean40 = list()
d_mat_mean50 = list()
d_mat_mean60 = list()
d_mat_mean90 = list()


error_analisis = '3_layer'


fig1 = figure(1)
ax = axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(d_mat_i, d_mat_i, '--', c='black', label='expected')
line_expected = Line2D([0], [0], color='black', ls='--')
for i in range(d_mat_i.size):
    if pDr[i] == -10:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='.', capsize=2, lw=1, c='purple')
        line10 = Line2D([0], [0], color='purple', ls='', marker='.')
        lmda10.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -20:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line20 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda20.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -30:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='^', capsize=2, lw=1, c='green')
        line30 = Line2D([0], [0], color='green', ls='', marker='^')
        lmda30.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -40:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='v', capsize=2, lw=1, c='orange')
        line40 = Line2D([0], [0], color='orange', ls='', marker='v')
        lmda40.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean40.append(d_mat_mean_i[i])
    elif pDr[i] == -50:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='1', capsize=2, lw=1, c='red')
        line50 = Line2D([0], [0], color='red', ls='', marker='1')
        lmda50.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -60:
            ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                        ls='', marker='1', capsize=2, lw=1, c='pink')
            line60 = Line2D([0], [0], color='pink', ls='', marker='2')
            lmda60.append(lmda[i])
            if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
                d_mat_mean60.append(d_mat_mean_i[i])
    elif pDr[i] == -90:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line90 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda90.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean90.append(d_mat_mean_i[i])
# for i in range(d_mat.size):
#     ax.annotate('(' + str(d_mat_pDr[i]) + ', ' + str(round(d_mat_mean[i], 1)) + ')',
#                 (d_mat[i], d_mat_mean[i])
#                 )
xlabel(r'$d_{sim}\ (\mu m)$')
ylabel(r'$d_{fit}\ (\mu m)$')
# xlim([d_mat[0], d_mat[-1]])
legend()  # loc='upper left')
# custom_lines = [line_expected, line10, line20, line30, line40, line50, line60]
# ax.legend(custom_lines, ['sim', -10, -20, -30, -40, -50, -60])
custom_lines = [line_expected, line60]
ax.legend(custom_lines, ['sim',-60])
savefig('./output/simulation_results/' + sample_dir + '/d_mat_fit_i_' + error_analisis + '.png')
# show()
# quit()
fig1 = figure(2)
ax = axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(d_mat_i, d_mat_i, '--', c='black', label='expected')
line_expected = Line2D([0], [0], color='black', ls='--')
for i in range(d_mat_i.size):
    if pDr[i] == -10:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='.', capsize=2, lw=1, c='purple')
        line10 = Line2D([0], [0], color='purple', ls='', marker='.')
        lmda10.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -20:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line20 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda20.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -30:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='^', capsize=2, lw=1, c='green')
        line30 = Line2D([0], [0], color='green', ls='', marker='^')
        lmda30.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -40:
        ax.errorbar(d_mat_m[i], d_mat_mean_m[i], yerr=d_mat_std_m[i],
                    ls='', marker='v', capsize=2, lw=1, c='orange')
        line40 = Line2D([0], [0], color='orange', ls='', marker='v')
        lmda40.append(lmda[i])
        if abs(d_mat_m[i] - d_mat_mean_m[i]) / d_mat_m[i] > 0.1:
            d_mat_mean40.append(d_mat_mean_m[i])
    elif pDr[i] == -50:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='1', capsize=2, lw=1, c='red')
        line50 = Line2D([0], [0], color='red', ls='', marker='1')
        lmda50.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -60:
            ax.errorbar(d_mat_m[i], d_mat_mean_m[i], yerr=d_mat_std_m[i],
                        ls='', marker='1', capsize=2, lw=1, c='pink')
            line60 = Line2D([0], [0], color='pink', ls='', marker='2')
            lmda60.append(lmda[i])
            if abs(d_mat_m[i] - d_mat_mean_m[i]) / d_mat_m[i] > 0.1:
                d_mat_mean60.append(d_mat_mean_m[i])
    elif pDr[i] == -90:
        ax.errorbar(d_mat_m[i], d_mat_mean_m[i], yerr=d_mat_std_m[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line90 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda90.append(lmda[i])
        if abs(d_mat_m[i] - d_mat_mean_m[i]) / d_mat_m[i] > 0.1:
            d_mat_mean90.append(d_mat_mean_m[i])
# for i in range(d_mat.size):
#     ax.annotate('(' + str(d_mat_pDr[i]) + ', ' + str(round(d_mat_mean[i], 1)) + ')',
#                 (d_mat[i], d_mat_mean[i])
#                 )
xlabel(r'$d_{sim}\ (\mu m)$')
ylabel(r'$d_{fit}\ (\mu m)$')
# xlim([d_mat[0], d_mat[-1]])
legend()  # loc='upper left')
# custom_lines = [line_expected, line10, line20, line30, line40, line50, line60]
# ax.legend(custom_lines, ['sim', -10, -20, -30, -40, -50, -60])
custom_lines = [line_expected, line60]
ax.legend(custom_lines, ['sim',-60])
savefig('./output/simulation_results/' + sample_dir + '/d_mat_fit_m_' + error_analisis + '.png')

fig1 = figure(3)
ax = axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(d_mat_i, d_mat_i, '--', c='black', label='expected')
line_expected = Line2D([0], [0], color='black', ls='--')
for i in range(d_mat_i.size):
    if pDr[i] == -10:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='.', capsize=2, lw=1, c='purple')
        line10 = Line2D([0], [0], color='purple', ls='', marker='.')
        lmda10.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -20:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line20 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda20.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -30:
        ax.errorbar(d_mat_i[i], d_mat_mean_i[i], yerr=d_mat_std_i[i],
                    ls='', marker='^', capsize=2, lw=1, c='green')
        line30 = Line2D([0], [0], color='green', ls='', marker='^')
        lmda30.append(lmda[i])
        if abs(d_mat_i[i] - d_mat_mean_i[i]) / d_mat_i[i] > 0.1:
            d_mat_mean10.append(d_mat_mean_i[i])
    elif pDr[i] == -40:
        ax.errorbar(d_mat_o[i], d_mat_mean_o[i], yerr=d_mat_std_o[i],
                    ls='', marker='v', capsize=2, lw=1, c='orange')
        line40 = Line2D([0], [0], color='orange', ls='', marker='v')
        lmda40.append(lmda[i])
        if abs(d_mat_o[i] - d_mat_mean_o[i]) / d_mat_o[i] > 0.1:
            d_mat_mean40.append(d_mat_mean_o[i])
    elif pDr[i] == -60:
            ax.errorbar(d_mat_o[i], d_mat_mean_o[i], yerr=d_mat_std_o[i],
                        ls='', marker='1', capsize=2, lw=1, c='pink')
            line60 = Line2D([0], [0], color='pink', ls='', marker='2')
            lmda60.append(lmda[i])
            if abs(d_mat_o[i] - d_mat_mean_o[i]) / d_mat_o[i] > 0.1:
                d_mat_mean60.append(d_mat_mean_o[i])
    elif pDr[i] == -90:
        ax.errorbar(d_mat_i[i], d_mat_mean_o[i], yerr=d_mat_std_o[i],
                    ls='', marker='+', capsize=2, lw=1, c='blue')
        line90 = Line2D([0], [0], color='blue', ls='', marker='+')
        lmda90.append(lmda[i])
        if abs(d_mat_o[i] - d_mat_mean_o[i]) / d_mat_o[i] > 0.1:
            d_mat_mean90.append(d_mat_mean_o[i])
# for i in range(d_mat.size):
#     ax.annotate('(' + str(d_mat_pDr[i]) + ', ' + str(round(d_mat_mean[i], 1)) + ')',
#                 (d_mat[i], d_mat_mean[i])
#                 )
xlabel(r'$d_{sim}\ (\mu m)$')
ylabel(r'$d_{fit}\ (\mu m)$')
# xlim([d_mat[0], d_mat[-1]])
legend()  # loc='upper left')
# custom_lines = [line_expected, line10, line20, line30, line40, line50, line60]
# ax.legend(custom_lines, ['sim', -10, -20, -30, -40, -50, -60])
custom_lines = [line_expected, line60]
ax.legend(custom_lines, ['sim', -60])
savefig('./output/simulation_results/' + sample_dir + '/d_mat_fit_o_' + error_analisis + '.png')

t_ref, E_ref = read_1file('./output/refs/-60_ref.txt')
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
n_i, k_i = nk_from_eps(mean(e_s_i), mean(e_inf_i), mean(tau_i), f_ref)
n_i_mean, k_i_mean = nk_from_eps(mean(e_s_mean_i), mean(e_inf_mean_i), mean(tau_mean_i), f_ref)
n_i_std, k_i_std = nk_from_eps(mean(e_s_std_i), mean(e_inf_std_i), mean(tau_std_i), f_ref)
n_m, k_m = nk_from_eps(mean(e_s_m), mean(e_inf_m), mean(tau_m), f_ref)
n_m_mean, k_m_mean = nk_from_eps(mean(e_s_mean_m), mean(e_inf_mean_m), mean(tau_mean_m), f_ref)
n_m_std, k_m_std = nk_from_eps(mean(e_s_std_m), mean(e_inf_std_m), mean(tau_std_m), f_ref)
n_o, k_o = nk_from_eps(mean(e_s_o), mean(e_inf_o), mean(tau_o), f_ref)
n_o_mean, k_o_mean = nk_from_eps(mean(e_s_mean_o), mean(e_inf_mean_o), mean(tau_mean_o), f_ref)
n_o_std, k_o_std = nk_from_eps(mean(e_s_std_o), mean(e_inf_std_o), mean(tau_std_o), f_ref)

f_min_idx, f_max_idx = f_min_max_idx(f_ref)
f_ref *= 1e-12

figure(11)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref, n_i, '-', c='black', label='expected')
ax.plot(f_ref, n_i_mean)
# ax.plot(f_ref[f_min_idx:f_max_idx], n_i[f_min_idx:f_max_idx], '-', c='black', label='expected')
# ax.plot(f_ref[f_min_idx:f_max_idx], n_i_mean[f_min_idx:f_max_idx])
# ax.errorbar(f_ref[f_min_idx:f_max_idx], n_i_mean[f_min_idx:f_max_idx], yerr=n_i_std[f_min_idx:f_max_idx])
figure(12)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref[f_min_idx:f_max_idx], k_i[f_min_idx:f_max_idx], '-', c='black', label='expected')
ax.errorbar(f_ref[f_min_idx:f_max_idx], k_i_mean[f_min_idx:f_max_idx], yerr=k_i_std[f_min_idx:f_max_idx])
figure(21)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref[f_min_idx:f_max_idx], n_m[f_min_idx:f_max_idx], '-', c='black', label='expected')
ax.errorbar(f_ref[f_min_idx:f_max_idx], n_m_mean[f_min_idx:f_max_idx], yerr=n_m_std[f_min_idx:f_max_idx])
figure(22)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref[f_min_idx:f_max_idx], k_m[f_min_idx:f_max_idx], '-', c='black', label='expected')
ax.errorbar(f_ref[f_min_idx:f_max_idx], k_m_mean[f_min_idx:f_max_idx], yerr=k_m_std[f_min_idx:f_max_idx])
figure(31)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref[f_min_idx:f_max_idx], n_o[f_min_idx:f_max_idx], '-', c='black', label='expected')
ax.errorbar(f_ref[f_min_idx:f_max_idx], n_o_mean[f_min_idx:f_max_idx], yerr=n_o_std[f_min_idx:f_max_idx])
figure(32)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(f_ref[f_min_idx:f_max_idx], k_o[f_min_idx:f_max_idx], '-', c='black', label='expected')
ax.errorbar(f_ref[f_min_idx:f_max_idx], k_o_mean[f_min_idx:f_max_idx], yerr=k_o_std[f_min_idx:f_max_idx])
show()


# figure(33)
# ax = axes()
# # ax.set_xscale('log')
# ax.set_yscale('log')
# medians_d_mat = array((median(d_mat_mean10), median(d_mat_mean20), median(d_mat_mean30), median(d_mat_mean40), median(d_mat_mean50)))
# medians_lmda = array((median(lmda10), median(lmda20), median(lmda30), median(lmda40), median(lmda50)))
# d_mat_vs_lmda = array((medians_d_mat, medians_lmda))
# plot(medians_lmda, medians_d_mat, '.')
# xlabel(r'$\lambda\ (\mu m)$')
# ylabel(r'$d_{lim}\ (\mu m)$')
#
# figure(44)
# ax = axes()
# # ax.set_xscale('log')
# ax.set_yscale('log')
# medians_d_mat = array((median(d_mat_mean10), median(d_mat_mean20), median(d_mat_mean30), median(d_mat_mean40), median(d_mat_mean50)))
# medians_lmda = array((median(lmda10), median(lmda20), median(lmda30), median(lmda40), median(lmda50)))
# d_mat_vs_lmda = array((medians_d_mat, medians_lmda))
# plot(medians_lmda, medians_d_mat, '.')
# xlabel(r'$\lambda\ (\mu m)$')
# ylabel(r'$d_{lim}\ (\mu m)$')
#
# # show()
# # quit()
#
#
# fig22 = figure(22)
# ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
# plot(d_mat_i, lmda, '.')
# # show()
# # quit()
#
#
# freqs = arange(100) * 1e10  # Hz
# # print(freqs)
# # quit()
# n_sim_i, k_sim_i = nk_from_eps(mean(e_s_i), mean(e_inf_i), mean(tau_i), freqs)
# n_fit_i, k_fit_i = nk_from_eps(mean(e_s_mean_i), mean(e_inf_mean_i), mean(tau_mean_i), freqs)
# n_fit_upp, k_fit_upp = nk_from_eps(mean(e_s_mean_i + e_s_std_i), mean(e_inf_mean_i + e_inf_std_i), mean(tau_mean_i + tau_std_i), freqs)
# n_fit_dwn, k_fit_dwn = nk_from_eps(mean(e_s_mean_i - e_s_std_i), mean(e_inf_mean_i - e_inf_std_i), mean(tau_mean_i - tau_std_i), freqs)
# freqs *= 1e-12  # THz
#
# fig2 = figure(2)
# ax = axes()
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# ax.plot(freqs, n_sim_i, 'r-', label='expected')
# ax.plot(freqs, n_fit_i, label='fitted')
# # ax.errorbar(freqs, n_fit, yerr=(n_fit_upp, n_fit_dwn), label='fitted')
# xlabel(r'$f\ (THz)$')
# # ylabel(r'$n$')
# xlim([freqs[0], freqs[-1]])
# legend(loc='upper left')
# savefig('./output/n_fit_' + error_analisis + '.png')
#
# fig3 = figure(3)
# ax = axes()
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# ax.plot(freqs, k_sim_i, 'r-', label='expected')
# ax.plot(freqs, k_fit_i, label='fitted')
# # ax.errorbar(freqs, k_fit, yerr=(k_fit_upp, k_fit_dwn), label='fitted')
# xlabel(r'$f\ (THz)$')
# # ylabel(r'$n$')
# xlim([freqs[0], freqs[-1]])
# legend(loc='upper left')
# savefig('./output/k_fit_' + error_analisis + '.png')
#
#
# fig4 = figure(4)
# ax = axes()
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# ax.plot(d_mat_std_i, f_cutoff, 'o')  # , label='expected')
# # ax.plot(d_mat, pDr, 'b.', label='fitted')
# xlabel(r'$\frac{\sigma_{d_{fit}}}{d_{sim}}$')
# ylabel(r'$f\ (THz)$')
# # xlim([(d_mat_std / d_mat)[0], (d_mat_std / d_mat)[-1]])
# # legend(loc='upper left')
#
# fig5 = figure(5)
# # ax = axes()
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# ax.plot(d_mat_std_i, pDr, 'o')  # , label='expected')
# # ax.plot(d_mat, pDr, 'b.', label='fitted')
# xlabel(r'$\frac{\sigma_{d_{fit}}}{d_{sim}}$')
# ylabel('PDR (dB)')
# # xlim([(d_mat_std / d_mat)[0], (d_mat_std / d_mat)[-1]])
# # legend(loc='upper left')
# show()
# quit()
#
# print('\\begin{table}[H]')
# print('\t\\centering\n\t\\begin{tabular}{rrrr}')
# print('\t\t$d_{sim}\ (\mu m)$ & pDr (dB) & $\epsilon_{s}\ (sim)$ & $\epsilon_{s}\ (fit)$ & $\\frac{\Delta\epsilon_{s}}{\epsilon_{s}}$ \\\\')
# for i in range(d_mat.size - 1):
#     print('\t\t' + str(round(d_mat[i], 1)), '&',  abs(d_mat_pDr[i]), '&', round(e_s[i], 3), '&', round(e_s_mean[i], 3), '$\pm$', round(e_s_std[i], 3), '&', round(abs(e_s[i] - e_s_mean[i])/e_s[i], 3), '\\\\')
# print('\t\t' + str(round(d_mat[-1], 1)), '&',  abs(d_mat_pDr[-1]), '&', round(e_s[-1], 3), '&', round(e_s_mean[-1], 3), '$\pm$', round(e_s_std[-1], 3), '&', round(abs(e_s[-1] - e_s_mean[-1])/e_s[-1], 3))
# print('\t\\end{tabular}\n\t\\caption{$\\epsilon_s$}\n\t\\label{tab:e_s_' + error_analisis + '}')
# print('\\end{table}')
# print()
#
# print('\\begin{table}[H]')
# print('\t\\centering\n\t\\begin{tabular}{rrrr}')
# print('\t\t$d_{sim}\ (\\mu m)$ & pDr (dB) & $\\epsilon_{\\infty}\ (sim)$ & $\epsilon_{\\infty}\ (fit)$ & $\\frac{\Delta\epsilon_{\inf}}{\epsilon_{\inf}}$ \\\\')
# for i in range(d_mat.size - 1):
#     print('\t\t' + str(round(d_mat[i], 1)), '&',  abs(d_mat_pDr[i]), '&', round(e_inf[i], 3), '&', round(e_inf_mean[i], 3), '$\pm$', round(e_inf_std[i], 3), '&', round(abs(e_inf[i] - e_inf_mean[i])/e_inf[i], 3), '\\\\')
# print('\t\t' + str(round(d_mat[-1], 1)), '&',  abs(d_mat_pDr[-1]), '&', round(e_inf[-1], 3), '&', round(e_inf_mean[-1], 3), '$\pm$', round(e_inf_std[-1], 3), '&', round(abs(e_inf[-1] - e_inf_mean[-1])/e_inf[-1], 3))
# print('\t\\end{tabular}\n\t\\caption{$\\epsilon_{\\infty}$}\n\t\\label{tab:e_inf_' + error_analisis + '}')
# print('\\end{table}')
# print()
#
# print('\\begin{table}[H]')
# print('\t\\centering\n\t\\begin{tabular}{rrrr}')
# print('\t\t$d_{sim}\ (\mu m)$ & pDr (dB) & $\\tau\ (fs) (sim)$ & $\\tau\ (fs) (fit)$ & $\\frac{\Delta\\tau}{\\tau}$ \\\\')
# for i in range(d_mat.size - 1):
#     print('\t\t' + str(round(d_mat[i], 1)), '&',  abs(d_mat_pDr[i]), '&', round(tau[i]*1e15, 3), '&', round(tau_mean[i]*1e15, 3), '$\pm$', round(tau_std[i]*1e15, 3), '&', round(abs(tau[i] - tau_mean[i])/tau[i], 3), '\\\\')
# print('\t\t' + str(round(d_mat[-1], 1)), '&',  abs(d_mat_pDr[-1]), '&', round(tau[-1]*1e15, 3), '&', round(tau_mean[-1]*1e15, 3), '$\pm$', round(tau_std[-1]*1e15, 3), '&', round(abs(tau[-1] - tau_mean[-1])/tau[-1], 3))
# print('\t\\end{tabular}\n\t\\caption{$\\tau$}\n\t\\label{tab:tau_' + error_analisis + '}')
# print('\\end{table}')
#
# print('\\begin{figure}[H]')
# print('\t\\centering')
# print('\t\\begin{subfigure}[t]{0.5\\textwidth}')
# print('\t\t\\centering')
# print('\t\t\\includegraphics[scale=0.35]{n_fit_' + error_analisis + '.png}')
# print('\t\t\\caption{Refractive index}')
# print('\t\t\\label{fig:n_fit_' + error_analisis + '}')
# print('\t\\end{subfigure}%')
# print('\t~')
# print('\t\\begin{subfigure}[t]{0.5\\textwidth}')
# print('\t\t\\centering')
# print('\t\t\\includegraphics[scale=0.35]{k_fit_' + error_analisis + '.png}')
# print('\t\t\\caption{Extinction coefficient}')
# print('\t\t\\label{fig:k_fit_' + error_analisis + '}')
# print('\t\\end{subfigure}')
# print('\t\\caption{Resulting refractive index for $\\epsilon_s$, $\\epsilon_{\\infty}$ and $\\tau$ used in the trace simulations}')
# print('\t\\label{fig:optical_params_' + error_analisis + '}')
# print('\\end{figure}')
#
# print('\\begin{figure}[H]')
# print('\t\\centering')
# print('\t\\includegraphics[scale=0.75]{d_mat_fit_' + error_analisis + '.png}')
# print('\t\\caption{Fitted values for thicknesses.}')
# print('\t\\label{fig:d_mat_fit_' + error_analisis + '}')
# print('\\end{figure}')
#
# show()
