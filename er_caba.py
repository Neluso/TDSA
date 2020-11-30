from TDSA import *


t_ref, E_ref = read_1file('./data/er_caba/ref_stitch.txt')
t_sam, E_sam = read_1file('./data/er_caba/sam_stitch_ajustado.txt')
t_ref *= 1e-12
t_sam *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

absH_w = abs(E_sam_w / E_ref_w)
angH_w = jepsen_unwrap(t_ref, E_ref, t_sam, E_sam)

figure(1)
plot(f_ref, prettyfy(E_ref_w, amax(abs(E_ref_w))), label='ref')
plot(f_sam, prettyfy(E_sam_w, amax(abs(E_ref_w))), label='sam')
legend()

figure(2)
plot(f_ref, absH_w)

figure(3)
plot(f_ref, angH_w)

show()
