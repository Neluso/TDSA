from DPS_functions import *
from read_data import read_data
from tkinter import *
from tkinter.filedialog import askopenfilenames


def spectra_plot():
    Tk().withdraw()
    filenames = askopenfilenames(initialdir='./data', title='Select reference data')
    t_s, E_s, is_error = read_data(filenames)
    return 0