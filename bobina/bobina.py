from TDSA import *


t_sam, E_sam = read_1file('./rollo_cafe_2.txt')
t_sam2, E_sam2 = read_1file('./rollo_cafe_1_lamina.txt')

delta_t_sam = mean(diff(t_sam))
enlargement = 10 * E_sam.size
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_sam * arange(1, enlargement + 1)))
E_sam2 = zero_padding(E_sam2, 0, enlargement)
t_sam2 = concatenate((t_sam2, t_sam2[-1] * ones(enlargement) + delta_t_sam * arange(1, enlargement + 1)))

t_sam *= 1e-12
t_sam2 *= 1e-12
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_sam2, E_sam_w2 = fourier_analysis(t_sam2, E_sam2)
plot(f_sam, toDb(E_sam_w), label='bobina')
plot(f_sam2, toDb(E_sam_w2), label='trozo')
xlabel(r'f (THz)')
ylabel('(dB)')
legend()
show()
