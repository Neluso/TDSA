from TDSA_main_package import *


step = 10  # GHz
max_freq = 15
step /= 1e3
freqs = arange(max_freq / step + 1) * step  # from 0 to 10 THz
times = arange(freqs.size) / (freqs.size * mean(diff(freqs)))  # ps

spectrum = freqs * exp(- 3 * freqs)
# spectrum += flip(arange(1, freqs.size + 1)) / sum(spectrum)
print(freqs[argmax(spectrum)])
spectrum /= sum(spectrum)
# plot(freqs, spectrum)
# show()
# quit()
phase = exp(1j * 2 * pi * random.rand())
c_spectr = list()
for i in range(spectrum.size):
    c_spectr.append(complex(spectrum[i]) * phase)

c_spectr = array(c_spectr)
c_spectr += (-1 + 2 * random.rand(c_spectr.size)) * 1e-4  # noise floor at 80 dB

roll_point = int(round(freqs.size / 5))

E_sim = roll(irfft(c_spectr, times.size), roll_point)
E_sim /= sum(abs(E_sim))

plot(times, E_sim)

show()
