from TDSA import *


# t_ref_st, E_ref_st = read_1file('./latiguillo_30/ref_stitch.txt')
# t_ref_1, E_ref_1 = read_1file('./latiguillo_30/ref_1.txt')
# t_sam_st, E_sam_st = read_1file('./latiguillo_30/sam_stitch.txt')
# t_sam_1, E_sam_1 = read_1file('./latiguillo_30/sam_1.txt')
#
# f_ref, E_ref_w = fourier_analysis(t_ref_st, E_ref_st)
# f_sam, E_sam_w = fourier_analysis(t_sam_st, E_sam_st)
#
# figure(1)
# plot(t_ref_1, E_ref_1)
# plot(t_sam_1, E_sam_1)
#
# figure(2)
# plot(t_ref_st, E_ref_st)
# plot(t_sam_st, E_sam_st)
#
#
# H_w = E_sam_w/E_ref_w
#
# figure(3)
# plot(f_ref, prettyfy(abs(E_ref_w), max(abs(E_ref_w))))
# plot(f_sam, prettyfy(abs(E_sam_w), max(abs(E_ref_w))))
#
# figure(4)
# plot(f_ref, abs(H_w))
#
# show()


t_ref, E_ref = read_1file('./data/test.txt')
figure(1)
plot(t_ref, E_ref)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
figure(2)
plot(f_ref, toDb_0(E_ref_w))
show()
