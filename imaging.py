# This script takes the data from imaging_data directory named as 'PosV_X_PosH_Y.txt', where X is the row position and
# Y the column position. Also a 'ref.txt' must be measured and left in the same directory.


from read_data import read_data
from DPS_functions import *
from numpy import zeros, sum, abs, array, linspace, transpose
from aux_functions import *
import os
from matplotlib.pyplot import show, imshow, figure, colorbar, xlabel, ylabel, close, savefig, plot, title
from tkinter.ttk import *
from tkinter import *
from tkinter import messagebox


def imaging(show_plots, hor_offset, ver_offset, resolution, temp_window):
    
    files = os.listdir('./imaging_data/')
    popup = Toplevel()
    popup.geometry('390x65')
    popup.title('TDSA Building image')
    Label(popup).grid(row=0, column=0)
    Label(popup).grid(row=1, column=0)
    progress_bar = Progressbar(popup, orient="horizontal", length=380, mode="determinate", maximum=len(files), value=0)
    progress_bar.grid(row=1, column=1)
    popup.pack_slaves()

    t_ref, E_ref, is_error = read_data('./imaging_data/ref.txt')  # support/matrix without T-Ink (THz-Ink: i.e. lactose)
    if is_error:
        return 0
    
    t_aux_ref = t_ref  # _aux_ref keeps copy of _ref during windowing for main for
    E_aux_ref = E_ref
    t_aux_j = t_ref  # _aux_j junk windowing recovering
    E_aux_j = E_ref
    
    if temp_window == 'blackman_harris':
        t_ref, E_ref, t_aux_j, E_aux_j = bh_windowing(t_ref, E_ref, t_aux_j, E_aux_j)
    if temp_window == 'chebisehev':
        t_ref, E_ref, t_aux_j, E_aux_j = cheb_windowing(t_ref, E_ref, t_aux_j, E_aux_j)
    if temp_window == 'hann':
        t_ref, E_ref, t_aux_j, E_aux_j = hann_windowing(t_ref, E_ref, t_aux_j, E_aux_j)
    if temp_window == 'force_exp':
        t_ref, E_ref, t_aux_j, E_aux_j = force_exp_windowing(t_ref, E_ref, t_aux_j, E_aux_j)
    if temp_window == 'tukey':
        t_ref, E_ref, t_aux_j, E_aux_j = tukey_windowing(t_aux_ref, E_aux_ref, t_aux_j, E_aux_j)
    
    t_ref *= 1e-12  # seconds
    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)
    row_max = 1
    col_max = 1
    offset = array([ver_offset, hor_offset])  # Vertical, Horizontal
    
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    E_ref_w = abs(E_ref_w)
    f_tol = 0.1e12  # Hz
    f_min = 0.53e12 - f_tol  # Hz possible lactose peaks at 0.53e12, and 1.376e12
    f_max = 0.53e12 + f_tol  # Hz
    f_min_idx, f_max_idx = f_range(f_ref, f_min, f_max)
    
    pixel_data = list()
    for file in files:
        popup.update()
        progress_bar['value'] += 1
        if file != 'ref.txt':
            t_sam, E_sam, is_error = read_data('./imaging_data/' + file)
            if is_error:
                messagebox.showerror('Error opening ' + file)
                continue
            if progress_bar['value'] == 1:
                if temp_window == 'blackman_harris':
                    t_aux_j, E_aux_j, t_sam, E_sam = bh_windowing(t_aux_ref, E_aux_ref, t_sam, E_sam)
                if temp_window == 'chebisehev':
                    t_aux_j, E_aux_j, t_sam, E_sam = cheb_windowing(t_aux_ref, E_aux_ref, t_sam, E_sam)
                if temp_window == 'hann':
                    t_aux_j, E_aux_j, t_sam, E_sam = hann_windowing(t_aux_ref, E_aux_ref, t_sam, E_sam)
                if temp_window == 'force_exp':
                    t_aux_j, E_aux_j, t_sam, E_sam = force_exp_windowing(t_aux_ref, E_aux_ref, t_sam, E_sam)
                if temp_window == 'tukey':
                    t_aux_j, E_aux_j, t_sam, E_sam = tukey_windowing(t_aux_ref, E_aux_ref, t_sam, E_sam)
            t_sam *= 1e-12  # seconds
            f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
            E_sam_w = abs(E_sam_w)
            H_w = E_sam_w / E_ref_w

            # Data ordering
            file = file.replace('.txt', '')
            row_pos = int(round
                          ((float(file.split('_')[1]) - offset[0]) / resolution))  # row position
            if row_pos > row_max:
                row_max = row_pos
            col_pos = int(round((float(file.split('_')[3]) - offset[1]) / resolution))  # column position
            if col_pos > col_max:
                col_max = col_pos
            alpha = sum(H_w[f_min_idx:f_max_idx])  # non-dimensional --- sum(abs(E_sam)) / sum(abs(E_ref))
            pixel_data.append((row_pos, col_pos, alpha))
    
    if row_max > 1:
        row_max += 1
    if col_max > 1:
        col_max += 1
    
    data = zeros((row_max, col_max))
    for items in pixel_data:
        row_pos = int(items[0])
        col_pos = int(items[1])
        alpha = float(items[2])
        data[row_pos, col_pos] = alpha
    
    row_length = row_max * resolution  # mm
    col_length = col_max * resolution  # mm
    figure(1)
    imshow(data,
           cmap='bone',  # 'gray'
           origin='lower',
           extent=(col_length + offset[1], offset[1], offset[0], row_length + offset[0])
           )
    xlabel('mm')
    ylabel('mm')
    colorbar()
    popup.destroy()
    
    if not os.path.exists('./output/'):
        os.mkdir('./output/')
    if not os.path.exists('./output/imaging/'):
        os.mkdir('./output/imaging/')
    savefig('./output/imaging/image.svg', format='svg')
    
    if show_plots:
        show()
    
    close('all')
    
    return 0
