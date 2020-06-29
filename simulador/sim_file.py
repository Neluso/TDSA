from TDS_constants import *
from numpy import *


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layer_number = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
n_subs = 1.25  # - 0.000 * 1j  # substrate refractive index -- cork 2.0
