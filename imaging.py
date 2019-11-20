# TODO not exactly, review...
# This script takes the data from imaging_data directory named as 'X_Y.txt', where X is the row position and Y the
# column position. Also a 'ref.txt' must be measured and left in the same directory.
# X and Y are just position indicators, spatial resolution must be added to the script in the data visualization.
# They must be labeled from 1 to N.


from read_data import read_data
from DPS_functions import *
from numpy import zeros, sum, abs, array, linspace, transpose
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure, colorbar, xlabel, ylabel, close, savefig, plot, title


def imaging(show_plots, hor_offset, ver_offset):
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
    offset = array([ver_offset, hor_offset])  # Vertical, Horizontal
    
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    E_ref_w = abs(E_ref_w)
    f_tol = 0.3e12  # Hz
    f_min = 1.376e12 - f_tol  # Hz possible peaks at 0.53e12, and 1.376e12
    f_max = 1.376e12 + f_tol  # Hz
    f_min_idx, f_max_idx = f_range(f_ref, f_min, f_max)
    E_ref_w_range = sum(E_ref_w[f_min_idx:f_max_idx])
    E_ref_range = sum(abs(E_ref))
    
    pixel_data = list()
    alpha_list = list()
    # print('Reading data')

    for file in files:
        if file != 'ref.txt':
            t_sam, E_sam, is_error = read_data('./imaging_data/' + file)
            if is_error:
                print('Error opening ' + file)
                return 0
            t_sam *= 1e-12  # seconds
            f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
            E_sam_w = abs(E_sam_w)
            E_sam_w_range = sum(E_sam_w[f_min_idx:f_max_idx])
            E_sam_range = sum(abs(E_sam))

            # Data ordering
            file = file.replace('.txt', '')
            row_pos = int((float(file.split('_')[1]) - offset[0]) / resolution)  # row position
            if row_pos > row_max:
                row_max = row_pos
            col_pos = int((float(file.split('_')[3]) - offset[1]) / resolution)  # column position
            if col_pos > col_max:
                col_max = col_pos
            alpha = E_sam_w_range / E_ref_w_range  # adimensional
            alpha_list.append(alpha)
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
    imshow(data, cmap='gray', origin='lower', extent=(col_length + offset[0] + offset[1], offset[1], offset[0], row_length + offset[0]))
    xlabel('mm')
    ylabel('mm')
    colorbar()
    savefig('./output/imaging/image.svg', format='svg')
    
    figure(2)
    data = transpose(array(data))[:][0]
    alpha_list = array(alpha_list)
    x_data = linspace(1, data.size, data.size)
    # plot(x_data, data, '.')
    p1 = polyfit(x_data, data, 1)
    y_data = p1[0] * x_data + p1[1]
    plot(x_data, y_data - data, '.')
    plot(x_data, zeros(x_data.size), '--')
    title('Dispersion')
    
    if not os.path.exists('./output/imaging/'):
        os.mkdir('./output/imaging/')
    savefig('./output/imaging/image.svg', format='svg')
    
    if show_plots:
        show()
    
    close('all')
    
    return 0
