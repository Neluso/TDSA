from numpy import *


def theta_2(n1, n2, theta1):
    return arcsin(n1 * sin(theta1) / n2)


def reflection(n1, n2, k1=0, k2=0, theta1=0):  # theta1 in radians
    theta2 = theta_2(n1, n2, theta1)
    n1_cplx = n1 - 1j * k1
    n2_cplx = n2 -1j * k2
    return (n2_cplx - n1_cplx) / (n1_cplx + n2_cplx)