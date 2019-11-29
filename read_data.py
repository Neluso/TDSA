from numpy import array
from tkinter import messagebox


def read_data(filenames):
    if filenames == None:
        messagebox.showerror('Error', 'No file selected')
        return 0, 0, True
    x_data = list()
    y_data = list()
    for filename in filenames:
        if filename == '':
            messagebox.showerror('Error', 'No file selected')
            return 0, 0, True
        try:
            fh = open(filename)
        except:
            messagebox.showerror('Error', 'Error opening ' + filename.split('/')[-1])
            continue
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
        x_data.append(array(x))
        y_data.append(array(y))
    return array(x_data), array(y_data), False
