# This script takes the data from imaging_data directory named as 'X_Y.txt', where X is the row position and Y the column
# position. Also a 'ref.txt' must be measured and left in the same directory.
# X and Y are just position indicators, spatial resolution must be added to the script in the data visualization.
# They must be labeled from 1 to N.


from read_data import read_data
from DPS_functions import *
from jepsen_index import jepsen_index
from numpy import argmax, real, zeros, sum, abs
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure, colorbar
from tqdm import tqdm


thick = 1.95e-3  # mm to m
files = os.listdir('./imaging_data/')
t_ref, E_ref = read_data('./imaging_data/ref.txt')
E_ref_sum = sum(abs(E_ref))
t_ref *= 1e-12  # seconds
nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)
row_max = 1
col_max = 1


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_min_idx, f_max_idx = f_min_max_idx(f_ref)
f_maximum_idx = argmax(E_ref_w)


pixel_data = list()
print('Reading data')
for file in tqdm(files):
    if file != 'ref.txt':
        t_sam, E_sam = read_data('./imaging_data/' + file)
        t_sam *= 1e-12  # seconds
        n, alpha_w = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness=thick)
        file = file.replace('.txt', '')
        row_pos = int(file.split('_')[0])  # row position
        if row_pos > row_max:
            row_max = row_pos
        col_pos = int(file.split('_')[1])  # column position
        if col_pos > col_max:
            col_max = col_pos
        pixel_data.append((row_pos, col_pos, sum(abs(E_sam)) / E_ref_sum))  # 0.01 * real(alpha_w[f_maximum_idx])))


print('Building image')
data = zeros((row_max, col_max))
for items in tqdm(pixel_data):
    row_pos = int(items[0]) - 1
    col_pos = int(items[1]) - 1
    alpha = 100 * items[2]
    data[row_pos, col_pos] = alpha
print('Image has been built')


print('Plotting...')
figure(1)
imshow(data)  # TODO add resolution to the data visualization
colorbar()
show()
