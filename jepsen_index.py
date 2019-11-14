# This script extracts the refractive index of a sample measured in a transmision THz-TDS using the algorithm proposed by
# Jepsen, P.U. J Infrared Milli Terahz Waves (2019) 40: 395. https://doi.org/10.1007/s10762-019-00578-0


from numpy import array, unwrap, pi, argmax, abs, exp, log10, angle, amax, polyfit, where, ones, log
from numpy.fft import rfft, rfftfreq
from DPS_functions import fourier_analysis
from TDS_constants import *


def toDb(x):
    return 20 * log10(abs(x))


def prettyfy(x, norm):
    return toDb(x / norm)


def nextpow2(i):
    n = 1
    while n < i:
        n *= 2
    return n


def centre_loc(E_data):  # finds the centre of the pulse based on
    t_0_pos = argmax(abs(E_data))
    return t_0_pos


def fourier_analysis(t_data, E_data, nSamp):
    samp_int = (t_data[1] - t_data[0])  # seconds
    E_data_w = rfft(E_data, n=nSamp)
    f_data = rfftfreq(nSamp, d=samp_int)  # Hz
    return f_data, E_data_w


def refractive_index(frq, delta_phi, thick):
    n = ones(delta_phi.size)
    for i in range(delta_phi.size):
        if frq[i] == 0:
            continue
        n[i] += c_0 * delta_phi[i] / (2 * pi * frq[i] * thick)
    return n


def n_quocient(ref_ind):
    n_quo = list()
    for i in range(ref_ind.size):
        n_quo.append(((ref_ind[i] + n_aire) ** 2) / (4 * ref_ind[i] * n_aire))
    return array(n_quo)


def alpha_w(ref_ind, H_0, thick):
    n_quo = n_quocient(ref_ind)
    alf = list()
    thick = - (2 / thick)  # m^-1
    for i in range(ref_ind.size):
        alf.append(log(H_0[i] * n_quo[i]))
    return thick * array(alf)


def jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness):  # TODO only return what is needed

    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)

    # Step 1: Finding the centre of the pulse to get t_0ref and t_0sam
    pos_t_0ref = centre_loc(E_ref)
    t_0ref = t_ref[pos_t_0ref]
    pos_t_0sam = centre_loc(E_ref)
    t_0sam = t_ref[pos_t_0sam]

    # Step 2: Fourier transform of measures
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
    H_w = list()  # complex transfer function initialisation
    nFFT = f_ref.size
    for i in range(nFFT):
        H_w.append(E_sam_w[i] / E_ref_w[i])
    H_w = array(H_w)

    # Step 3: Calculate reduced phases
    phi_0_ref = list()
    phi_0_sam = list()
    for i in range(nFFT):
        phi_0_ref.append(2 * pi * f_ref[i] * t_0ref)
        phi_0_sam.append(2 * pi * f_sam[i] * t_0sam)
    phi_0_ref = array(phi_0_ref)
    phi_0_sam = array(phi_0_sam)
    phi_0_ref_red = list()
    phi_0_sam_red = list()
    for i in range(nFFT):
        phi_0_ref_red.append(E_ref_w[i] * exp(- 1j * phi_0_ref[i]))
        phi_0_sam_red.append(E_sam_w[i] * exp(- 1j * phi_0_sam[i]))
    phi_0_ref_red = array(phi_0_ref_red)
    phi_0_sam_red = array(phi_0_sam_red)
    phi_0_ref_red = angle(phi_0_ref_red)
    phi_0_sam_red = angle(phi_0_sam_red)

    # Step 4: Unwrap the reduced phase difference
    delta_phi_0_red = abs(unwrap(phi_0_sam_red - phi_0_ref_red))

    # Step 5: Fit the unwrapped phase to a linear function
    f_min_idx = 1
    f_max_idx = 100
    f_min = f_ref[f_min_idx]
    f_max = f_ref[f_max_idx]
    for frq in f_ref:
        if frq <= 1e11:  # 0.1 THz
            f_min = frq
        if frq <= 1e12:  # 1 THz
            f_max = frq
    f_min_idx = where(f_ref == f_min)[0][0]
    f_max_idx = where(f_ref == f_max)[0][0]
    coefs = polyfit(f_ref[f_min_idx:f_max_idx], delta_phi_0_red[f_min_idx:f_max_idx], 2)

    delta_phi_0 = delta_phi_0_red - 2 * pi * ones(delta_phi_0_red.size) * int(coefs[2] / (2 * pi))
    delta_phi = abs(delta_phi_0 + (phi_0_sam - phi_0_ref))

    # Step 6.1: Obtaining the refractive index
    n = refractive_index(f_ref, delta_phi, thickness)

    # Step 6.2: Obtaining the absorption coefficient
    alpha_f = alpha_w(n, H_w, thickness)

    return f_ref, f_min_idx, f_max_idx, n, alpha_f
