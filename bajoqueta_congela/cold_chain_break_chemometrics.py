from TDSA import *
from scipy.stats import t, expon, norm
from scipy.optimize import curve_fit


def self_unwrap(t_sam, E_sam):
    t_sam_0 = centroid_E2(t_sam, E_sam)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    phi_sam_0 = 2 * pi * f_sam * t_sam_0
    phi_sam_0_red = angle(E_sam_w * exp(-1j * phi_sam_0))
    return unwrap(phi_sam_0_red)


def data_fit(x, mu, sigma):
    return norm.pdf(x, mu, sigma)


sigmas = list()
descong = list()

for i in range(8):
    # if i+1 == 5 or i+1 == 6: continue
    thick = open('./cold_chain_break/' + str(i + 1) + 'a_bajoqueta_rotura/thickness.txt')
    thick = float(thick.read())
    # thick = 1e3  # mm
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

        delta_t_ref = mean(diff(t_ref))
        enlargement = 9 * E_ref.size
        E_ref = zero_padding(E_ref, 0, enlargement)
        t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
        E_sam = zero_padding(E_sam, 0, enlargement)
        t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))

        f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
        
        # f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.15, 0.85)
        f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.25, 0.35)

        H_w = E_sam_w / E_ref_w
        H_w = H_w[f_min_idx:f_max_idx]
        f_ref = f_ref[f_min_idx:f_max_idx]
        E_ref_w = E_ref_w[f_min_idx:f_max_idx]
        filt_ref = abs(E_ref_w)/max(abs(E_ref_w))
        # plot(f_ref * 1e-12, filt_ref)
        # show()
        # quit()
        
        gH_w = gradient(H_w)
        aH_w = abs(H_w)
        agH_w = abs(gH_w)
        gaH_w = gradient(aH_w)
        ggH_w = gradient(gH_w)
        ggaH_w = gradient(gaH_w)
        gagH_w = gradient(abs(gH_w))
        aggH_w = abs(ggH_w)

        phi_unwrapped = jepsen_unwrap(t_ref, E_ref, t_sam, E_sam)
        gPhi = gradient(phi_unwrapped)
        agPhi = abs(gPhi)
        ggPhi = gradient(gPhi)
        gagPhi = gradient(agPhi)
        aggPhi = abs(ggPhi)
        
        hist_data = ggPhi[f_min_idx:f_max_idx]  # gPhi
        # hist_data = aggH_w
        
        bins_numb = int(hist_data.size / 2)
        
        hist_data_obj = histogram(hist_data, bins=bins_numb)

        sigmas.append(round(std(hist_data), 4))
        descong.append(j + 1)
        
        # figure(1993)
        # if i == 0:
        #     colr = 'blue'
        # elif i ==1:
        #     colr = 'orange'
        # elif i ==2:
        #     colr = ''
        # plot(j + 1, round(std(hist_data), 4), '.', label=str(i + 1), color=colr)

        figure(10 + i)
        hist(hist_data + 1*j, bins=bins_numb, label=r'$\sigma =$' + str(round(std(hist_data), 4)))
        legend()
        savefig('./cold_chain_break/output/hist_' + str(i + 1) + '.png')
        close()

        figure(i + 1)
        title('Mostra ' + str(i + 1))
        plot(f_ref, hist_data + 1*j, label=r'$\sigma =$' + str(round(std(hist_data), 4)))  # round(popt[1], 4)))
        legend()
        savefig('./cold_chain_break/output/phase_' + str(i + 1) + '.png')
        close()


descong = array(descong)
sigmas = array(sigmas)
# sigmas = toDb(sigmas)
p2 = polyfit(descong, sigmas, 2)
p1 = polyfit(descong, sigmas, 1)
figure(1993)
plot(descong, sigmas, '.')
plot(descong, p1[0] * descong + p1[1])
plot(arange(1, 5), p2[0] * arange(1, 5)**2 + p2[1] * arange(1, 5) + p2[2])
# legend()
show()
