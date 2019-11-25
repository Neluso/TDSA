# This script extracts the refractive index of a sample measured in a transmision THz-TDS using the algorithm proposed by
# Jepsen, P.U. J Infrared Milli Terahz Waves (2019) 40: 395. https://doi.org/10.1007/s10762-019-00578-0


from numpy import array, unwrap, pi, abs, exp, angle, polyfit, ones, log
from aux_functions import *
from DPS_functions import *
from TDS_constants import *


def refractive_index(frq, delta_phi, thick):
    n = ones(delta_phi.size)
    for i in range(delta_phi.size):
        if frq[i] == 0:
            continue
        n[i] += c_0 * delta_phi[i] / (2 * pi * frq[i] * thick)
    return n


def n_quocient(ref_ind):
    return ((ref_ind + n_aire) ** 2) / (4 * ref_ind * n_aire)


def alpha_w(ref_ind, H_0, thick):
    return - (2 / thick) * log(H_0 * n_quocient(ref_ind))  # m^-1


def jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness):  # Returns refractive index 'n' and absortion coeficient 'alpha_r'

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
    H_w = E_sam_w / E_ref_w  # complex transfer function

    # Step 3: Calculate reduced phases
    phi_0_ref = 2 * pi * f_ref * t_0ref
    phi_0_sam = 2 * pi * f_sam * t_0sam
    phi_0_ref_red = E_ref_w * exp(- 1j * phi_0_ref)
    phi_0_sam_red = E_sam_w * exp(- 1j * phi_0_sam)
    phi_0_ref_red = angle(phi_0_ref_red)
    phi_0_sam_red = angle(phi_0_sam_red)

    # Step 4: Unwrap the reduced phase difference
    delta_phi_0_red = abs(unwrap(phi_0_sam_red - phi_0_ref_red))

    # Step 5: Fit the unwrapped phase to a linear function
    f_min_idx, f_max_idx = f_min_max_idx(f_ref)

    coefs = polyfit(f_ref[f_min_idx:f_max_idx], delta_phi_0_red[f_min_idx:f_max_idx], 2)

    delta_phi_0 = delta_phi_0_red - 2 * pi * ones(delta_phi_0_red.size) * int(coefs[2] / (2 * pi))
    delta_phi = abs(delta_phi_0 + (phi_0_sam - phi_0_ref))

    # Step 6.1: Obtaining the refractive index
    n = refractive_index(f_ref, delta_phi, thickness)

    # Step 6.2: Obtaining the absorption coefficient in m^-1
    alpha_f = alpha_w(n, H_w, thickness)

    return n, alpha_f
