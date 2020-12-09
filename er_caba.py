from TDSA import *

for i in range(3):
    t_ref, E_ref = read_1file('./data/er_caba/v2/ref' + str(i+1) + '_1.txt')
    t_sam, E_sam = read_1file('./data/er_caba/v2/sam' + str(i+1) + '_1.txt')
    
    delta_t_ref = mean(diff(t_ref))
    enlargement = 10 * E_ref.size
    E_ref = zero_padding(E_ref, 0, enlargement)
    t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    E_sam = zero_padding(E_sam, 0, enlargement)
    t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
    
    
    # plot(t_ref, E_ref)
    # plot(t_sam, E_sam)
    # show()
    # quit()
    t_ref *= 1e-12
    t_sam *= 1e-12
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
    
    absH_w = abs(E_sam_w / E_ref_w)
    angH_w = jepsen_unwrap(t_ref, E_ref, t_sam, E_sam)
    
    f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.1, 0.4)
    f_ref *= 1e-12
    f_sam *= 1e-12
    
    figure(1)
    # if i == 0:
    #     plot(f_ref[f_min_idx:f_max_idx], prettyfy(E_ref_w, amax(abs(E_ref_w)))[f_min_idx:f_max_idx], label='ref')
    plot(f_sam, prettyfy(E_sam_w, amax(abs(E_ref_w))), label='sam')
    xlabel('f (THz)')
    # legend()
    
    figure(2)
    plot(f_ref[f_min_idx:f_max_idx], absH_w[f_min_idx:f_max_idx])
    ylim([0, 1.5])
    xlabel('f (THz)')
    
    figure(3)
    plot(f_ref[f_min_idx:f_max_idx], angH_w[f_min_idx:f_max_idx])
    xlabel('f (THz)')

show()
