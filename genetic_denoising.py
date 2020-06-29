from TDSA import *
from scipy import interpolate
from scipy.signal import savgol_filter


t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)


# Matlab runs idxs differently
def smooth(M, span):
    M_aux = M
    p_max = M.size
    for p in range(p_max):
        if p - span < 0:
            M_aux[p] = sum(M[:2 * p + 1]) / (2 * p + 1)
        elif span < p - 1 < p_max - 1 - span:
            M_aux[p] = sum(M[p - span:p + span]) / (2 * span + 1)
        elif p + span > p_max - 1:
            M_aux[p] = sum(M[2 * p - p_max - 1:p_max - 1]) / (2 * p_max - 2 * p - 1)
    return M_aux


def score(M, E_sam, E_ref_w):
    return std(E_sam - irfft(M*E_ref_w))


def genetic_deno(t_ref, E_ref, t_sam, E_sam):

    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

    num_pop = 200
    num_pop_cross = int(num_pop/2)
    num_pop_new = num_pop - num_pop_cross
    max_iter = 500
    noise_level_cutoff = 2e-5  # noise level for iteration cutoff
    smoothing_span = 4
    smoothing_order = 1
    p_mut = 0.2  # % of individuals which will suffer mutation
    mut_sev = 0.15  # 'severity' of mutation
    M_i = dict()  # Guesses
    S_i = dict()  # Scores

    M_raw = E_sam_w / E_ref_w

    # M_smooth = smooth(M_raw, smoothing_span)
    M_smooth = savgol_filter(abs(M_raw), 2*smoothing_span+1, smoothing_order)
    M_smooth = M_smooth * exp(1j * savgol_filter(unwrap(angle(M_raw)), 2*smoothing_span+1, smoothing_order))
    # plot(fftshift(f_ref), fftshift(abs(M_raw)), 'o')
    # plot(fftshift(f_ref), fftshift(abs(M_smooth)), 'o')

    # figure(1)
    # plot(f_ref, abs(M_raw), label='raw', lw=0.5)
    # plot(f_sam, abs(M_smooth), label='smooth', lw=0.5)
    # figure(2)
    # plot(f_sam, unwrap(angle(M_raw)), label='raw', lw=0.5)
    # plot(f_sam, unwrap(angle(M_smooth)), label='smooth', lw=0.5)
    # show()
    # quit()


    # Initilization and 1st evaluation
    for i in range(num_pop):
        idxs = [0, M_raw.size-1]
        for idx in random.randint(0, M_raw.size, int(0.5*M_raw.size)):
            if idx not in idxs:
                idxs.append(idx)
        idxs = sort(idxs)
        rand_abs = (abs(M_raw[idxs]) - abs(M_smooth[idxs])) * random.rand(len(idxs)) + abs(M_smooth[idxs])
        # TODO reducir la fase con Jepsen
        # pos_t_0ref = centre_loc(E_ref)
        # t_0ref = t_ref[pos_t_0ref]
        # pos_t_0sam = centre_loc(E_ref)
        # t_0sam = t_ref[pos_t_0sam]
        # phi_0_ref = 2 * pi * f_ref * t_0ref
        # phi_0_sam = 2 * pi * f_sam * t_0sam
        # phi_0_ref_red = E_ref_w * exp(- 1j * phi_0_ref)
        # phi_0_sam_red = E_sam_w * exp(- 1j * phi_0_sam)
        # phi_0_ref_red = angle(phi_0_ref_red)
        # phi_0_sam_red = angle(phi_0_sam_red)
        # delta_phi_0_red = unwrap(phi_0_sam_red - phi_0_ref_red)
        ##############
        rand_arg = (unwrap(angle(M_raw[idxs])) - unwrap(angle(M_smooth[idxs]))) * random.rand(len(idxs)) + unwrap(angle(M_smooth[idxs]))
        # rand_arg = (delta_phi_0_red[idxs] - unwrap(angle(M_smooth[idxs]))) * random.rand(len(idxs)) + unwrap(angle(M_smooth[idxs]))
        abs_M_itpl = interpolate.interp1d(f_ref[idxs], rand_abs, fill_value='extrapolate')
        angle_M_itpl = interpolate.interp1d(f_ref[idxs], rand_arg, fill_value='extrapolate')
        m_i = abs_M_itpl(f_ref) * exp(1j * angle_M_itpl(f_ref))

        M_i[i] = m_i
        S_i[i] = score(m_i, E_sam, E_ref_w)


    for i in range(max_iter):

        # Crossover
        key_list = list()
        for item in sorted(S_i.items(), key=lambda x: x[1]):
            key_list.append(item[0])

        if S_i[key_list[0]] <= noise_level_cutoff:
            break

        M_i_cross = dict()
        S_i_cross = dict()
        new_key = 0
        for key in key_list[:num_pop_cross]:
            M_i_cross[new_key] = M_i[key]
            S_i_cross[new_key] = S_i[key]
            new_key += 1

        M_i = M_i_cross
        S_i = S_i_cross

        key_list = list()
        for item in sorted(S_i.items(), key=lambda x: x[1]):
            key_list.append(item[0])

        # Mutation
        for j in range(-1, num_pop_cross - 1):
            m_i = 0.5*(abs(M_i[key_list[j]]) + abs(M_i[key_list[j + 1]])) * exp(1j * 0.5 * (unwrap(angle(M_i[key_list[j]])) + unwrap(angle(M_i[key_list[j + 1]]))))
            p_rand = random.rand()
            if p_rand <= p_mut:
                f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.1, 1.5)
                # mut_idx = random.randint(0, M_raw.size)
                mut_idx = random.randint(f_min_idx, f_max_idx)
                if mut_idx <= smoothing_span:
                    if mut_idx == 0:
                        idxs = [0]
                        idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
                    else:
                        idxs = [0, mut_idx]
                        idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
                elif mut_idx >= M_raw.size - smoothing_span:
                    idxs = list(range(mut_idx - smoothing_span))
                    idxs = idxs + [mut_idx, M_raw.size - 1]
                else:
                    idxs = list(range(mut_idx - smoothing_span))
                    idxs.append(mut_idx)
                    idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
                m_i[mut_idx] = m_i[mut_idx] * (2 * random.rand() - 1) * mut_sev
                abs_M_itpl = interpolate.interp1d(f_ref[idxs], abs(m_i[idxs]), fill_value='extrapolate')
                angle_M_itpl = interpolate.interp1d(f_ref[idxs], unwrap(angle(m_i[idxs])), fill_value='extrapolate')
                m_i = abs_M_itpl(f_ref) * exp(1j * angle_M_itpl(f_ref))

            M_i[new_key] = m_i
            S_i[new_key] = score(m_i, E_sam, E_ref_w)
            new_key += 1

        key_list = list()
        for item in sorted(S_i.items(), key=lambda x: x[1]):
            key_list.append(item[0])


        # figure(1)
        # plot(f_sam, abs(M_fit), label='fit_'+str(i), lw=0.5)
        # figure(2)
        # plot(f_sam, unwrap(angle(M_fit)), label='fit_'+str(i), lw=0.5)

        # Mutation
        # for j in range(num_pop):
        #     p_rand = random.rand()
        #     if p_rand <= p_mut:
        #         f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.1, 1.0)
        #         mut_idx = random.randint(0, M_raw.size)
        #         if mut_idx <= smoothing_span:
        #             if mut_idx == 0:
        #                 idxs = [0]
        #                 idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
        #             else:
        #                 idxs = [0, mut_idx]
        #                 idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
        #         elif mut_idx >= M_raw.size - smoothing_span:
        #             idxs = list(range(mut_idx - smoothing_span))
        #             idxs = idxs + [mut_idx, M_raw.size-1]
        #         else:
        #             idxs = list(range(mut_idx - smoothing_span))
        #             idxs.append(mut_idx)
        #             idxs = idxs + list(range(mut_idx + smoothing_span, M_raw.size))
        #         m_j = M_i[j]
        #         m_j[mut_idx] = m_j[mut_idx] * (2 * random.rand() - 1) * mut_sev
        #         abs_M_itpl = interpolate.interp1d(f_ref[idxs], abs(m_j[idxs]), fill_value='extrapolate')
        #         angle_M_itpl = interpolate.interp1d(f_ref[idxs], unwrap(angle(m_j[idxs])), fill_value='extrapolate')
        #         m_i = abs_M_itpl(f_ref) * exp(1j * angle_M_itpl(f_ref))
        #         M_i[j] = m_i
        #         S_i[j] = score(m_i, E_sam, E_ref_w)

        key_list = list()
        for item in sorted(S_i.items(), key=lambda x: x[1]):
            key_list.append(item[0])
        M_fit = M_i[key_list[0]]
        print('Iteration ', i + 1, 'of', max_iter, '--- best score =', score(M_fit, E_sam, E_ref_w))
        # figure(20)
        # plot(i + 1, score(M_fit, E_sam, E_ref_w), 'bo', lw=0.5)


    key_list = list()
    for item in sorted(S_i.items(), key=lambda x: x[1]):
        key_list.append(item[0])
    M_fit = M_i[key_list[0]]

    return irfft(M_fit * E_ref_w)

    # figure()
    # # plot(t_ref, E_ref)
    # plot(t_sam, E_sam, label='raw')
    # plot(t_sam, real(irfft(M_smooth * E_ref_w)), label='smooth')
    # plot(t_sam, real(irfft(M_fit * E_ref_w)), label='fit')
    #
    # legend()
    #
    #
    # figure(1)
    # plot(f_sam, abs(M_fit), label='fit', lw=0.5)
    #
    # xlim([0, 1])
    # legend()
    #
    # figure(2)
    # plot(f_sam, unwrap(angle(M_fit)), label='fit', lw=0.5)
    #
    # xlim([0, 1])
    # legend()
    #
    # show()


# genetic_deno(t_ref, E_ref, t_sam, E_sam)