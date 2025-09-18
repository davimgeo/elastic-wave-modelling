import numpy as np
import matplotlib.pyplot as plt

SAVE_GEOM = False

SRC_DEPTH = 10
REC_DEPTH = 6

nb       = 100
nx, nz   = 1150, 648
nxx, nzz = nx + 2 * nb, nz + 2 * nb

salt = np.fromfile(
    "data/input/model_vp_2d_1150x648.bin", dtype=np.float32, count=nx*nz
).reshape([nz, nx], order='F')

# rec_pos / rec_depth
rec = np.zeros((nx, 2))

rec[:, 0] = np.arange(0, nx)
rec[:, 1] = REC_DEPTH

# src_pos / src_depth
src = np.zeros((nx, 2))

src[:, 0] = nx // 2
src[:, 1] = SRC_DEPTH

plt.plot(rec[:,0], rec[:,1], 'b*')
plt.plot(src[:,0], src[:,1], 'ro')

plt.imshow(salt)

plt.show()

if SAVE_GEOM:
    np.savetxt(
     "data/input/receivers.txt", 
     rec,
     fmt='%.2f',
     delimiter=',', 
     header='rec_pos, rec_depth'
    )

    np.savetxt(
     "data/input/sources.txt", 
     src,
     fmt='%.2f',
     delimiter=',', 
     header='src_pos, src_depth'
    )



