import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

PATH = "data/output/snapshots/"

nb = 30
nx, nz = 1150, 648
nxx, nzz = nx + 2 * nb, nz + 2 * nb

salt = np.fromfile("data/output/vp.bin", dtype=np.float32, count=nxx * nzz).reshape([nzz, nxx], order='F')

snap = np.fromfile(
        "data/output/snapshots/vp_1210x708_tid_830.bin", dtype=np.float32, count=nxx*nzz
).reshape([nzz, nxx], order='F')

scale = 0.4 * np.std(snap)

plt.imshow(snap, cmap='Greys', vmin=-scale, vmax=scale)
plt.imshow(salt, alpha=0.5, cmap='gray')

plt.show()




