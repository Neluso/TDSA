# This script takes the data from imaging_data directory named as 'X_Y.txt', where X is the row position and Y the column
# position. Also a 'ref.txt' must be measured and left in the same directory.
# X and Y are just position indicators, spatial resolution must be added to the script in the data visualization.
# They must be labeled from 1 to N.


from read_data import read_data
from DPS_functions import *
from jepsen_index import jepsen_index
from numpy import argmax, real, ones
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure
from tqdm import tqdm


thick = 1.95e-3  # mm to m
files = os.listdir('./imaging_data/')
t_ref, E_ref = read_data('./imaging_data/ref.txt')
t_ref *= 1e-12  # seconds
nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_min_idx, f_max_idx = f_min_max_idx(f_ref)
f_maximum_idx = argmax(E_ref_w)


pixel_data = list()
for file in files:
    if file != 'ref.txt':
        t_sam, E_sam = read_data('./imaging_data/' + file)
        t_sam *= 1e-12  # seconds
        n, alpha_w = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness=thick)
        file = file.replace('.txt', '')
        row_pos = file.split('_')[0]  # row position
        col_pos = file.split('_')[1]  # column position
        pixel_data.append((row_pos, col_pos, 0.01 * real(alpha_w[f_maximum_idx])))


print('Building image')
# Square image
data = ones((int(len(pixel_data) / 2), int(len(pixel_data) / 2)))
for items in tqdm(pixel_data):
    row_pos = int(items[0]) - 1
    col_pos = int(items[1]) - 1
    alpha = items[2]
    data[row_pos, col_pos] = alpha
print('Image has been built')


print('Plotting...')
figure(1)
imshow(data)  # TODO add resolution to the data visualization
show()
