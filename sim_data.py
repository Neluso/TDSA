from TDSA import *


step = 10  # GHz
max_freq = 15
step /= 1e3
freqs = arange(max_freq / step + 1) * step  # from 0 to 10 THz
times = arange(freqs.size) / (freqs.size * mean(diff(freqs)))  # ps

spectrum = freqs * exp(- 3 * freqs)
spectrum /= sum(spectrum)

phase = exp(1j * 2 * pi * random.rand())
c_spectr = list()
for i in range(spectrum.size):
    c_spectr.append(complex(spectrum[i]) * phase)

c_spectr = array(c_spectr)
c_spectr += (-1 + 2 * random.rand(c_spectr.size)) * fromDb(-80)  # noise floor at -60 dB

roll_point = int(round(freqs.size / 5))

E_sim = roll(irfft(c_spectr, times.size), roll_point)
E_sim /= sum(abs(E_sim))

sim_pulse_idx = centre_loc(E_sim)
pulse_dist = 0.50
material_pulse = pulse_dist * roll(E_sim, sim_pulse_idx - 100)
material_pulse += (1 - pulse_dist) * roll(material_pulse, 10)
plot(times, E_sim, lw=1)
plot(times, material_pulse, lw=1)

show()
