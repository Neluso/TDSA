from TDSA import *
import csv


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


d_mat = array(d_mat) * 1e6
d_mat_mean = array(d_mat_mean) * 1e6
d_mat_std = array(d_mat_std) * 1e6
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


fig1 = figure(1)
ax = axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(d_mat, d_mat, 'r--', label='expected')
ax.errorbar(d_mat, d_mat_mean, yerr=d_mat_std, label='fitted')
xlabel(r'$d_{sim}\ (\mu m)$')
ylabel(r'$d_{fit}\ (\mu m)$')
xlim([d_mat[0], d_mat[-1]])
legend(loc='upper left')

print('$d_{sim}\ (\mu m)$ & $\epsilon_{s}\ (sim)$ & $\epsilon_{s}\ (fit)$ & $\\frac{\Delta\epsilon_{s}}{\epsilon_{s}}$ \\\\')
for i in range(d_mat.size):
    print(round(d_mat[i], 1), '&', round(e_s[i], 3), '&', round(e_s_mean[i], 3), '$\pm$', round(e_s_std[i], 3), '&', round(abs(e_s[i] - e_s_mean[i])/e_s[i], 3), '\\\\')

print()
print('$d_{sim}\ (\mu m)$ & $\epsilon_{\inf}\ (sim)$ & $\epsilon_{\inf}\ (fit)$ & $\\frac{\Delta\epsilon_{\inf}}{\epsilon_{\inf}}$ \\\\')
for i in range(d_mat.size):
    print(round(d_mat[i], 1), '&', round(e_inf[i], 3), '&', round(e_inf_mean[i], 3), '$\pm$', round(e_inf_std[i], 3), '&', round(abs(e_inf[i] - e_inf_mean[i])/e_inf[i], 3), '\\\\')

print()
print('$d_{sim}\ (\mu m)$ & $\\tau\ (THz) (sim)$ & $\\tau\ (THz) (fit)$ & $\\frac{\Delta\\tau}{\\tau}$ \\\\')
for i in range(d_mat.size):
    print(round(d_mat[i], 1), '&', round(tau[i]*1e12, 3), '&', round(tau_mean[i]*1e12, 3), '$\pm$', round(tau_std[i]*1e12, 3), '&', round(abs(tau[i] - tau_mean[i])/tau[i], 3), '\\\\')


show()
