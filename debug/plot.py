import matplotlib.pyplot as plt
import numpy as np

PATH = "data/output/wavelet.txt"

ricker = np.loadtxt(PATH)

t = np.linspace(0, 1001*4e-3, len(ricker))

plt.plot(t, ricker)

plt.show()