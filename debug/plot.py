import matplotlib.pyplot as plt
import numpy as np

PATH = "data/output/snapshots/txx_501x501_tid_800.bin"
#PATH = "data/output/txx.bin"
#PATH = "data/input/model_rho_2d_1150x648.bin"
#PATH = "data/output/vp.bin"

nx, nz = 501, 501
#nx, nz = 1150, 648

model = np.fromfile(PATH, dtype=np.float32, count = nx * nz).reshape([nz, nx], order='F') 

plt.imshow(model, cmap='gray')
# plt.show()
plt.savefig("model.png")  # Save instead of showing
print("Plot saved as model.png")
