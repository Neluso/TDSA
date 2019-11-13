from numpy import unwrap, pi, linspace
from matplotlib.pyplot import plot, show


n = 10
phase = linspace(0, pi, num=n)
x = linspace(0, n, num=n)
phase[5:] += 2 * pi
plot(x, phase)
plot(x, unwrap(unwrap(phase)))
show()
