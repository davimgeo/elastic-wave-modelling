import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# PATH = "data/output/vp.bin"
#
# nb = 30
# nx, nz = 1150, 648
# nxx, nzz = nx + 2*nb, nz + 2*nb
#
# model = np.fromfile(PATH, dtype=np.float32, count=nxx * nzz).reshape([nzz, nxx], order='F')
#
# plt.imshow(model, cmap='Greys')
#
# plt.show()

PATH = "data/output/snapshots/"

nx, nz = 1150, 648

salt = np.fromfile(
        "data/input/model_vp_2d_1150x648.bin", dtype=np.float32, count=nx * nz
).reshape([nz, nx], order='F')

snap = np.fromfile(
        "data/output/snapshots/vp_1150x648_tid_1660.bin", dtype=np.float32, count=nx*nz
).reshape([nz, nx], order='F')

scale = 0.6 * np.std(snap)

plt.imshow(snap, cmap='Greys', vmin=-scale, vmax=scale)
plt.imshow(salt, alpha=0.5, cmap='gray')

plt.show()




