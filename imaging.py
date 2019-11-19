# This script takes the data from imaging_data directory named as 'X_Y.txt', where X is the row position and Y the column
# position. Also a 'ref.txt' must be measured and left in the same directory.
# X and Y are just position indicators, spatial resolution must be added to the script in the data visualization.
# They must be labeled from 1 to N.


from read_data import read_data
from DPS_functions import *
from numpy import zeros, sum, abs
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure, colorbar, xlabel, ylabel
from tqdm import tqdm


thick = 1.95e-3  # mm to m
files = os.listdir('./imaging_data/')
t_ref, E_ref = read_data('./imaging_data/ref.txt')  # support/matrix without T-Ink (THz-Ink: i.e. lactose)
t_ref *= 1e-12  # seconds
nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)
resolution = 1  # mm
row_max = 1
col_max = 1


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
E_ref_w = abs(E_ref_w)
f_min = 1.2e12  # Hz
f_max = 1.5e12  # Hz
f_min_idx, f_max_idx = f_range(f_ref, f_min, f_max)
E_ref_w_range = sum(E_ref_w[f_min_idx:f_max_idx])


pixel_data = list()
print('Reading data')
for file in tqdm(files):
    if file != 'ref.txt':
        t_sam, E_sam = read_data('./imaging_data/' + file)
        t_sam *= 1e-12  # seconds
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
        E_sam_w = abs(E_sam_w)
        E_sam_w_range = sum(E_sam_w[f_min_idx:f_max_idx])
        
        # Data ordering
        file = file.replace('.txt', '')
        row_pos = int(file.split('_')[0])  # row position
        if row_pos > row_max:
            row_max = row_pos
        col_pos = int(file.split('_')[1])  # column position
        if col_pos > col_max:
            col_max = col_pos
        alpha = E_sam_w_range / E_ref_w_range  # adimensional
        pixel_data.append((row_pos, col_pos, alpha))
print('Data has been read')


print('Building image')
data = zeros((row_max, col_max))
for items in tqdm(pixel_data):
    row_pos = int(items[0]) - 1
    col_pos = int(items[1]) - 1
    alpha = 100 * float(items[2])
    data[row_pos, col_pos] = alpha
print('Image has been built')


print('Plotting')
row_length = row_max * resolution  # mm
col_length = col_max * resolution  # mm
figure(1)
imshow(data, extent=(0, col_length, row_length, 0))
xlabel('mm')
ylabel('mm')
colorbar()
show()
