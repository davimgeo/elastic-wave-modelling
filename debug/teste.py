import numpy as np
import matplotlib.pyplot as plt

nb = 30
nx, nz = 1150, 648
nxx, nzz = nx + 2*nb, nz + 2*nb

salt = np.fromfile(
    "data/output/text.bin", dtype=np.float32, count=nxx*nzz
).reshape([nzz, nxx], order='F')

plt.imshow(salt)

plt.show()


