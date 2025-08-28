import matplotlib.pyplot as plt
import numpy as np

# PATH = "data/output/wavelet.txt"
#
# ricker = np.loadtxt(PATH)
#
# t = np.linspace(0, 1001*4e-3, len(ricker))
#
# plt.plot(t, ricker)
#
# plt.show() 
#
# PATH = "data/output/snapshots/vx_501x501_tid_495.bin"
PATH = "data/output/txx.bin"
nx, ny = 501, 501

model = np.fromfile(PATH, dtype=np.float32)
model = model.reshape((nx, ny), order='F') 

plt.imshow(model.T, cmap='gray') 
plt.show()
