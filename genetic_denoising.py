from TDSA import *
from scipy import interpolate


t_ref, E_ref = read_1file('./data/demo_data/test_ref.txt')
t_sam, E_sam = read_1file('./data/demo_data/test_sam.txt')
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)


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


def score(M):
    return std(E_sam - ifft(M*E_ref_w))


num_pop = 1000
num_pop_cross = int(num_pop/2)
num_pop_new = num_pop - num_pop_cross
max_iter = 1000
smoothing_span = 2
p_mut = 0.1  # % of individuals which will suffer mutation
M_i = dict()


M_raw = E_sam_w / E_ref_w
M_smooth = smooth(M_raw, smoothing_span)


# Initilization and 1st evaluation
for i in range(num_pop):
    idxs = [0, M_raw.size-1]
    for idx in random.randint(0, M_raw.size, int(0.5*M_raw.size)):
        if idx not in idxs:
            idxs.append(idx)
    idxs = sort(idxs)
    abs_M_itpl = interpolate.interp1d(f_ref[idxs], abs(M_smooth[idxs]), fill_value='extrapolate')
    # plot(fftshift(f_ref), fftshift(abs(M_raw)))
    # plot(fftshift(f_ref), fftshift(abs_M_itpl(f_ref)))
    # show()
    angle_M_itpl = interpolate.interp1d(f_ref[idxs], unwrap(angle(M_smooth[idxs])), fill_value='extrapolate')
    m_i = abs_M_itpl(f_ref) * exp(-1j * angle_M_itpl(f_ref))

    M_i[score(m_i)] = m_i


for i in trange(max_iter):

    # Crossover
    key_list = sorted(M_i.keys())
    key_list = key_list[:num_pop_cross]
    M_i_cross = dict()
    for key in key_list:
        M_i_cross[key] = M_i[key]
    M_i = M_i_cross

    for j in range(num_pop_cross - 1):
        m_i = 0.5*(M_i[key_list[j]] + M_i[key_list[j + 1]])
        M_i[score(m_i)] = m_i

    # Mutation


key_list = sorted(M_i.keys())
M_fit = M_i[key_list[0]]


figure()
plot(t_sam, E_sam, label='sam')
plot(t_sam, real(ifft(smooth(M_fit, smoothing_span) * E_ref_w)), label='fit')

legend()


figure()
plot(fftshift(f_ref), fftshift(abs(M_raw)), label='raw')
plot(fftshift(f_sam), fftshift(abs(smooth(M_fit, 5))), label='smooth')

xlim([0, 1])
legend()

show()
