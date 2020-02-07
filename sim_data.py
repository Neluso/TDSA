from TDSA import *


step = 10  # GHz
max_freq = 15
step /= 1e3
freqs = arange(max_freq / step + 1) * step  # from 0 to 10 THz
times = arange(freqs.size) / (freqs.size * mean(diff(freqs)))  # ps
spectrum = freqs * exp(- 3 * freqs)
t_ref, E_ref = read_1file('./data/sim_resources/transmision_ref.txt')
t_ref2, E_ref2 = read_1file('./data/sim_resources/refletion_ref.txt')
# t_ref, E_ref = read_1file('./data/sim_resources/mod_transmision_ref.txt')
# t_ref2, E_ref2 = read_1file('./data/sim_resources/mod_refletion_ref.txt')

plot(t_ref, E_ref)
plot(t_ref2, E_ref2)



show()
