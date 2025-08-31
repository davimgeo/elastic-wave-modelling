import os
import numpy as np
import matplotlib.pyplot as plt

PATH = "data/output/snapshots/"

nx, nz = 1150, 648

salt = np.fromfile(
        "data/input/salt_model/model_vp_2d_1150x648.bin", dtype=np.float32, count=nx * nz
).reshape([nz, nx], order='F')

snap = np.fromfile(
        "data/output/snapshots/txx_1150x648_tid_1494.bin", dtype=np.float32, count=nx*nz
).reshape([nz, nx], order='F')

scale = 0.8 * np.std(snap)

plt.imshow(snap, cmap='Greys', vmin=-scale, vmax=scale)
plt.imshow(salt, alpha=0.5)

plt.show()




