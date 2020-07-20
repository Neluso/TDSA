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

d_mat = list()
d_mat_mean = list()
d_mat_std = list()

e_s = list()
e_s_mean = list()
e_s_std = list()

e_inf = list()
e_inf_mean = list()
e_inf_std = list()

tau = list()
tau_mean = list()
tau_std = list()

d_air = list()
d_air_mean = list()
d_air_std = list()

f_cuttof = list()
pDr = list()

# rows = list()

with open('./output/resolution_limit.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        d_mat.append(float(row[0]))
        d_mat_mean.append(float(row[1]))
        d_mat_std.append(float(row[2]))
        e_s.append(float(row[3]))
        e_s_mean.append(float(row[4]))
        e_s_std.append(float(row[5]))
        e_inf.append(float(row[6]))
        e_inf_mean.append(float(row[7]))
        e_inf_std.append(float(row[8]))
        tau.append(float(row[9]))
        tau_mean.append(float(row[10]))
        tau_std.append(float(row[11]))
        d_air.append(float(row[12]))
        d_air_mean.append(float(row[13]))
        d_air_std.append(float(row[14]))
        f_cuttof.append(float(row[15]))
        pDr.append(float(row[16]))



d_mat = array(d_mat)  # * 1e6
d_mat_mean = array(d_mat_mean)  # * 1e6
d_mat_std = array(d_mat_std)  # * 1e6
d_mat_pDr = array(pDr)
e_s = array(e_s)
e_s_mean = array(e_s_mean)
e_s_std = array(e_s_std)
e_inf = array(e_inf)
e_inf_mean = array(e_inf_mean)
e_inf_std = array(e_inf_std)
tau = array(tau)
tau_mean = array(tau_mean)
tau_std = array(tau_std)
d_air = array(tau)
d_air_mean = array(tau_mean)
d_air_std = array(tau_std)
f_cuttof = array(f_cuttof)
pDr = array(pDr)


error_analisis = 'p100_100'


fig1 = figure(1)
ax = axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(d_mat, d_mat, 'r--', label='expected')
ax.errorbar(d_mat, d_mat_mean,
            yerr=d_mat_std, label='fitted',
            ls='', marker='.', capsize=2, lw=1
            # uplims=True, lolims=True
            )
# for i in range(d_mat.size):
#     ax.annotate('(' + str(d_mat_pDr[i]) + ', ' + str(round(d_mat_mean[i], 1)) + ')',
#                 (d_mat[i], d_mat_mean[i])
#                 )
xlabel(r'$d_{sim}\ (\mu m)$')
ylabel(r'$d_{fit}\ (\mu m)$')
# xlim([d_mat[0], d_mat[-1]])
legend(loc='upper left')
savefig('./output/d_mat_fit_' + error_analisis + '.png')
show()
quit()


freqs = arange(100) * 1e10
# print(freqs)
# quit()
n_sim, k_sim = nk_from_eps(mean(e_s), mean(e_inf), mean(tau), freqs)
n_fit, k_fit = nk_from_eps(mean(e_s_mean), mean(e_inf_mean), mean(tau_mean), freqs)
n_fit_upp, k_fit_upp = nk_from_eps(mean(e_s_mean + e_s_std), mean(e_inf_mean + e_inf_std), mean(tau_mean + tau_std), freqs)
n_fit_dwn, k_fit_dwn = nk_from_eps(mean(e_s_mean - e_s_std), mean(e_inf_mean - e_inf_std), mean(tau_mean - tau_std), freqs)
freqs *= 1e-12

fig2 = figure(2)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(freqs, n_sim, 'r-', label='expected')
ax.plot(freqs, n_fit, label='fitted')
# ax.errorbar(freqs, n_fit, yerr=(n_fit_upp, n_fit_dwn), label='fitted')
xlabel(r'$f\ (THz)$')
# ylabel(r'$n$')
xlim([freqs[0], freqs[-1]])
legend(loc='upper left')
savefig('./output/n_fit_' + error_analisis + '.png')

fig3 = figure(3)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(freqs, k_sim, 'r-', label='expected')
ax.plot(freqs, k_fit, label='fitted')
# ax.errorbar(freqs, k_fit, yerr=(k_fit_upp, k_fit_dwn), label='fitted')
xlabel(r'$f\ (THz)$')
# ylabel(r'$n$')
xlim([freqs[0], freqs[-1]])
legend(loc='upper left')
savefig('./output/k_fit_' + error_analisis + '.png')


fig4 = figure(4)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(d_mat_std, f_cuttof, 'o')  # , label='expected')
# ax.plot(d_mat, pDr, 'b.', label='fitted')
xlabel(r'$\frac{\sigma_{d_{fit}}}{d_{sim}}$')
ylabel(r'$f\ (THz)$')
# xlim([(d_mat_std / d_mat)[0], (d_mat_std / d_mat)[-1]])
# legend(loc='upper left')

fig5 = figure(5)
ax = axes()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.plot(d_mat_std, toDb(pDr), 'o')  # , label='expected')
# ax.plot(d_mat, pDr, 'b.', label='fitted')
xlabel(r'$\frac{\sigma_{d_{fit}}}{d_{sim}}$')
ylabel('PDR (dB)')
# xlim([(d_mat_std / d_mat)[0], (d_mat_std / d_mat)[-1]])
# legend(loc='upper left')

print('\\begin{table}[H]')
print('\t\\centering\n\t\\begin{tabular}{rrrr}')
print('\t\t$d_{sim}\ (\mu m)$ & $\epsilon_{s}\ (sim)$ & $\epsilon_{s}\ (fit)$ & $\\frac{\Delta\epsilon_{s}}{\epsilon_{s}}$ \\\\')
for i in range(d_mat.size - 1):
    print('\t\t' + str(round(d_mat[i], 1)), '&', round(e_s[i], 3), '&', round(e_s_mean[i], 3), '$\pm$', round(e_s_std[i], 3), '&', round(abs(e_s[i] - e_s_mean[i])/e_s[i], 3), '\\\\')
print('\t\t' + str(round(d_mat[-1], 1)), '&', round(e_s[-1], 3), '&', round(e_s_mean[-1], 3), '$\pm$', round(e_s_std[-1], 3), '&', round(abs(e_s[-1] - e_s_mean[-1])/e_s[-1], 3))
print('\t\\end{tabular}\n\t\\caption{$\\epsilon_s$}\n\t\\label{tab:e_s_' + error_analisis + '}')
print('\\end{table}')
print()

print('\\begin{table}[H]')
print('\t\\centering\n\t\\begin{tabular}{rrrr}')
print('\t\t$d_{sim}\ (\\mu m)$ & $\\epsilon_{\\infty}\ (sim)$ & $\epsilon_{\\infty}\ (fit)$ & $\\frac{\Delta\epsilon_{\inf}}{\epsilon_{\inf}}$ \\\\')
for i in range(d_mat.size - 1):
    print('\t\t' + str(round(d_mat[i], 1)), '&', round(e_inf[i], 3), '&', round(e_inf_mean[i], 3), '$\pm$', round(e_inf_std[i], 3), '&', round(abs(e_inf[i] - e_inf_mean[i])/e_inf[i], 3), '\\\\')
print('\t\t' + str(round(d_mat[-1], 1)), '&', round(e_inf[-1], 3), '&', round(e_inf_mean[-1], 3), '$\pm$', round(e_inf_std[-1], 3), '&', round(abs(e_inf[-1] - e_inf_mean[-1])/e_inf[-1], 3))
print('\t\\end{tabular}\n\t\\caption{$\\epsilon_{\\infty}$}\n\t\\label{tab:e_inf_' + error_analisis + '}')
print('\\end{table}')
print()

print('\\begin{table}[H]')
print('\t\\centering\n\t\\begin{tabular}{rrrr}')
print('\t\t$d_{sim}\ (\mu m)$ & $\\tau\ (fs) (sim)$ & $\\tau\ (fs) (fit)$ & $\\frac{\Delta\\tau}{\\tau}$ \\\\')
for i in range(d_mat.size - 1):
    print('\t\t' + str(round(d_mat[i], 1)), '&', round(tau[i]*1e15, 3), '&', round(tau_mean[i]*1e15, 3), '$\pm$', round(tau_std[i]*1e15, 3), '&', round(abs(tau[i] - tau_mean[i])/tau[i], 3), '\\\\')
print('\t\t' + str(round(d_mat[-1], 1)), '&', round(tau[-1]*1e15, 3), '&', round(tau_mean[-1]*1e15, 3), '$\pm$', round(tau_std[-1]*1e15, 3), '&', round(abs(tau[-1] - tau_mean[-1])/tau[-1], 3))
print('\t\\end{tabular}\n\t\\caption{$\\tau$}\n\t\\label{tab:tau_' + error_analisis + '}')
print('\\end{table}')

print('\\begin{figure}[H]')
print('\t\\centering')
print('\t\\begin{subfigure}[t]{0.5\\textwidth}')
print('\t\t\\centering')
print('\t\t\\includegraphics[scale=0.35]{n_fit_' + error_analisis + '.png}')
print('\t\t\\caption{Refractive index}')
print('\t\t\\label{fig:n_fit_' + error_analisis + '}')
print('\t\\end{subfigure}%')
print('\t~')
print('\t\\begin{subfigure}[t]{0.5\\textwidth}')
print('\t\t\\centering')
print('\t\t\\includegraphics[scale=0.35]{k_fit_' + error_analisis + '.png}')
print('\t\t\\caption{Extinction coefficient}')
print('\t\t\\label{fig:k_fit_' + error_analisis + '}')
print('\t\\end{subfigure}')
print('\t\\caption{Resulting refractive index for $\\epsilon_s$, $\\epsilon_{\\infty}$ and $\\tau$ used in the trace simulations}')
print('\t\\label{fig:optical_params_' + error_analisis + '}')
print('\\end{figure}')

print('\\begin{figure}[H]')
print('\t\\centering')
print('\t\\includegraphics[scale=0.75]{d_mat_fit_' + error_analisis + '.png}')
print('\t\\caption{Fitted values for thicknesses.}')
print('\t\\label{fig:d_mat_fit_' + error_analisis + '}')
print('\\end{figure}')

show()
