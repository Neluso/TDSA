from numpy import array
from tkinter import messagebox


# TODO improve data reading: read more sample data to perform statistical analysis


def read_data(filename):
    if filename == '':
        messagebox.showerror('Error', 'No file selected')
        return 0, 0, True
    try:
        fh = open(filename)
    except:
        messagebox.showerror('Error', 'Error opening ' + filename.split('/')[-1])
        return 0, 0, True
    data = fh.read()
    data = data.split('\n')
    x = list()
    y = list()
    for item in data:
        item = item.split(',')
        if item == ['']:
            break
        x.append(float(item[0]))
        y.append(float(item[1]))
    x = array(x)
    y = array(y)
    return x, y, False
