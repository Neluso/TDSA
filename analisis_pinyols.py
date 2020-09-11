from TDSA import *
from genetic_denoising import genetic_deno


for i in range(7):
    muestra_j = str(i + 1)
    dh = open('./pinyol_olives/data/muestra_' + str(muestra_j) + '/datos.txt')
    # print(dh.read().split('\n')[0].split(' '))
    # quit()
    medible = True
    try:
        meh1, thick, meh2, delta_thick, meh3 = dh.read().split('\n')[0].split(' ')
        thick = float(thick) * 1e-6
        delta_thick = float(delta_thick) * 1e-6
        print(i + 1)
    except:
        print('No medible', i + 1)
        medible = False
    
    n_idx, alpha_r = list(), list()
    
    for i in range(3):
        t_ref, E_ref = read_1file('./pinyol_olives/data/muestra_' + str(muestra_j) + '/ref_' + str(i + 1) + '.txt')
        t_sam, E_sam = read_1file('./pinyol_olives/data/muestra_' + str(muestra_j) + '/sam_' + str(i + 1) + '.txt')
        t_ref *= 1e-12
        t_sam *= 1e-12
        if medible:
            n_idx_aux, alpha_r_aux, n_idx_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, thick)
            n_idx.append(n_idx_aux)
            alpha_r.append(alpha_r_aux)
        
        f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
        # E_sam_w = genetic_deno(t_ref, E_ref, t_sam, E_sam)
        
        figure(1)
        if i == 0:
            plot(t_ref * 1e12, E_ref)
        plot(t_sam * 1e12, E_sam)
        savefig('./pinyol_olives/data/muestra_' + str(muestra_j) + '/trazas.jpg', format='JPG')
        clf()
        
        f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.1, 1.2)
        f_ref = f_ref[f_min_idx:f_max_idx]
        f_sam = f_sam[f_min_idx:f_max_idx]
        E_ref_w = E_ref_w[f_min_idx:f_max_idx]
        E_sam_w = E_sam_w[f_min_idx:f_max_idx]
        
        figure(2)
        if i == 0:
            plot(f_ref * 1e-12, prettyfy(E_ref_w, max(abs(E_ref_w))))
        plot(f_sam * 1e-12, prettyfy(E_sam_w, max(abs(E_ref_w))))
        savefig('./pinyol_olives/data/muestra_' + str(muestra_j) + '/espectros.jpg', format='JPG')
        clf()
    
    if medible:
        n_idx = array(n_idx)
        alpha_r = array(alpha_r)
        n_idx_mean = mean(n_idx, axis=0)
        alpha_r_mean = mean(alpha_r, axis=0)
        n_idx_std = std(n_idx, axis=0)
        alpha_r_std = std(alpha_r, axis=0)
        n_idx_mean = n_idx_mean[f_min_idx:f_max_idx]
        alpha_r_mean = alpha_r_mean[f_min_idx:f_max_idx]
        n_idx_std = n_idx_std[f_min_idx:f_max_idx]
        alpha_r_std = alpha_r_std[f_min_idx:f_max_idx]
        figure(3)
        errorbar(f_ref * 1e-12, n_idx_mean, yerr=n_idx_std)
        savefig('./pinyol_olives/data/muestra_' + str(muestra_j) + '/indice.jpg', format='JPG')
        clf()
        figure(4)
        errorbar(f_ref * 1e-12, alpha_r_mean * 1e-2, yerr=alpha_r_std * 1e-2)
        savefig('./pinyol_olives/data/muestra_' + str(muestra_j) + '/alpha.jpg', format='JPG')
        clf()
