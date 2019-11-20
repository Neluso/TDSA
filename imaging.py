# This script takes the data from imaging_data directory named as 'X_Y.txt', where X is the row position and Y the column
# position. Also a 'ref.txt' must be measured and left in the same directory.
# X and Y are just position indicators, spatial resolution must be added to the script in the data visualization.
# They must be labeled from 1 to N.


from read_data import read_data
from DPS_functions import *
from numpy import zeros, sum, abs, array
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure, colorbar, xlabel, ylabel, close, savefig, plot, axvline


def imaging(show_plots):
    files = os.listdir('./imaging_data/')
    t_ref, E_ref, is_error = read_data('./imaging_data/ref.txt')  # support/matrix without T-Ink (THz-Ink: i.e. lactose)
    if is_error:
        return 0
    t_ref *= 1e-12  # seconds
    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)
    resolution = 0.5  # mm
    row_max = 1
    col_max = 1
    offset = array([30, 20])  # Vertical, Horizontal
    
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    E_ref_w = abs(E_ref_w)
    f_tol = 0.5e12  # Hz
    f_min = 1.376e12 - f_tol  # 1.376e12 + f_tol  # Hz possible peaks at 0.53e12, and 1.376e12
    f_max = 1.376e12 + f_tol  # 1.376e12 - f_tol  # Hz
    f_min_idx, f_max_idx = f_range(f_ref, f_min, f_max)
    E_ref_w_range = sum(E_ref_w[f_min_idx:f_max_idx])
    figure(2)
    plot(f_ref[f_min_idx:f_max_idx], toDb(E_ref_w[f_min_idx:f_max_idx]), lw=7)
    axvline(1.376e12 + f_tol, c='r', ls='--', lw=0.5)
    axvline(1.376e12 - f_tol, c='r', ls='--', lw=0.5)
    
    pixel_data = list()
    # print('Reading data')
    for file in files:
        if file != 'ref.txt':
            t_sam, E_sam, is_error = read_data('./imaging_data/' + file)
            if is_error:
                return 0
            t_sam *= 1e-12  # seconds
            f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
            E_sam_w = abs(E_sam_w)
            E_sam_w_range = sum(E_sam_w[f_min_idx:f_max_idx])
            
            plot(f_sam[f_min_idx:f_max_idx], toDb(E_sam_w[f_min_idx:f_max_idx]), linewidth=0.5)
            
            # Data ordering
            file = file.replace('.txt', '')
            row_pos = int(2 * (float(file.split('_')[1]) - offset[0]))  # row position
            if row_pos > row_max:
                row_max = row_pos
            col_pos = int(2 * (float(file.split('_')[3]) - offset[1]))  # column position
            if col_pos > col_max:
                col_max = col_pos
            alpha = E_sam_w_range / E_ref_w_range  # adimensional
            pixel_data.append((row_pos, col_pos, alpha))
    # print('Data has been read')
    
    # print('Building image')
    data = zeros((row_max, col_max))
    for items in pixel_data:
        row_pos = int(items[0]) - 1
        col_pos = int(items[1]) - 1
        alpha = float(items[2])  # * 100
        data[row_pos, col_pos] = alpha
    # print('Image has been built')
    
    # print('Plotting')
    row_length = row_max * resolution  # mm
    col_length = col_max * resolution  # mm
    figure(1)
    imshow(data, origin='lower', extent=(offset[0] + offset[1], offset[1], offset[0], row_length + offset[0]))
    xlabel('mm')
    ylabel('mm')
    colorbar()
    
    if not os.path.exists('./output/imaging/'):
        os.mkdir('./output/imaging/')
    savefig('./output/imaging/image.svg', format='svg')
    
    if show_plots:
        show()
    
    close('all')
    
    return 0
