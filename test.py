from TDSA import *


t_ref, E_ref = read_1file('./data/for_testing/PP100/ref2.txt')
t_sam, E_sam = read_1file('./data/for_testing/PP100/sam2.txt')
# t_ref, E_ref = read_1file('./data/for_testing/PP70_PET30/ref3.txt')
# t_sam, E_sam = read_1file('./data/for_testing/PP70_PET30/sam3.txt')
t_ref *= 1e-12
t_sam *= 1e-12
# t_ref, E_ref = read_1file('./data/demo_data/ref.txt')
# t_sam, E_sam = read_1file('./data/demo_data/sam.txt')
# t_ref, E_ref = read_1file('./data/for_testing/PP70_PET30/ref4.txt')
# t_sam, E_sam = read_1file('./data/for_testing/1_1_PP.txt')
# plot(t_sam, E_sam)
# show()
# quit()


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

figure(1)
plot(f_ref, toDb(E_ref_w))
plot(f_sam, toDb(E_sam_w))

figure(2)
plot(f_ref, abs(E_sam_w / E_ref_w))
xlabel(r'$f\ (THz)$')
ylabel(r'$|H(\omega)|$')
xlim([0, 1.6])
ylim([0, 1.6])

figure(3)
n, alpha, n_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, 2.02e-3)
plot(f_ref, n)
figure(4)
plot(f_ref, alpha)
show()
