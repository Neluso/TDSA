from tkinter import Tk
from tkinter.filedialog import askopenfilename


Tk().withdraw()
filename = askopenfilename(initialdir='./data',title='Selecciona fichero de referencia')
print(filename.split('/')[-1])