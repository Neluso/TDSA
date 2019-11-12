# This script extracts the refractive index of a sample measured in a transmision THz-TDS using the algorithm proposed by
# Jepsen, P.U. J Infrared Milli Terahz Waves (2019) 40: 395. https://doi.org/10.1007/s10762-019-00578-0


from numpy import array, unwrap, pi, argmax, abs, exp, log10, angle, amax
from numpy.fft import rfft, rfftfreq


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
    samp_int = (t_data[1] - t_data[0])*1e-12  # seconds
    E_data_w = rfft(E_data, n=nSamp)
    f_data = rfftfreq(nSamp, d=samp_int)  # Hz
    return f_data, E_data_w


def jepsen_index(t_ref, E_ref, t_sam, E_sam):  # returns

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
    H_w = list()  # transfer function
    nFFT = f_ref.size
    for i in range(nFFT):
        H_w.append(E_ref_w[i] - E_sam_w[i])
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
    delta_phi_0_red = unwrap(phi_0_sam_red - phi_0_ref_red)

    return f_ref, prettyfy(E_ref_w, amax(E_ref_w)), f_sam, prettyfy(E_sam_w, amax(E_ref_w)), H_w / amax(H_w), delta_phi_0_red
