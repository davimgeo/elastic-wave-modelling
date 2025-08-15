import matplotlib.pyplot as plt
import numpy as np

PATH = "data/output/wavelet.txt"

ricker = np.loadtxt(PATH)

t = np.linspace(0, 1001*4e-3, len(ricker))

plt.plot(t, ricker)

plt.show() 

#PATH = "data/input/model_rho_2d_1150x648.bin"
#nx, ny = 1150, 648 

#model = np.fromfile(PATH, dtype=np.float32)
#model = model.reshape((nx, ny), order='F')  # reshape in Fortran order

#plt.imshow(model.T, cmap='viridis') 
#plt.show()
