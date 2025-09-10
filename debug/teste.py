import numpy as np
import matplotlib.pyplot as plt

SAVE_GEOM = True

SRC_DEPTH = 10
REC_DEPTH = 6

nb       = 100
nx, nz   = 1150, 648
nxx, nzz = nx + 2 * nb, nz + 2 * nb

salt = np.fromfile(
    "data/input/model_vp_2d_1150x648.bin", dtype=np.float32, count=nx*nz
).reshape([nz, nx], order='F')

# columns:
# rec_pos / rec_depth / src_pos / src_depth 
geom = np.zeros((nx, 4))

geom[:, 0] = np.arange(0, nx)
geom[:, 1] = REC_DEPTH
geom[:, 2] = nx // 2
geom[:, 3] = SRC_DEPTH

plt.plot(geom[:,0], geom[:,1], 'b*')
plt.plot(geom[:,2], geom[:,3], 'ro')

plt.imshow(salt)

plt.show()

if SAVE_GEOM:
    np.savetxt(
        "data/input/receivers.txt", 
        geom,
        fmt='%.2f',
        delimiter=',', 
        header='rec_pos, rec_depth, src_pos, src_depth'
    )



