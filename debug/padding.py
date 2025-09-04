import numpy as np
import matplotlib.pyplot as plt

nx, nz = 501, 501

nb = 30

nxx = nx + 2 * nb
nzz = nz + 2 * nb

z0  = 200
amp = 50
k   = (2 * np.pi) / 100

model = np.zeros((nzz, nxx))

x = np.arange(nx)    

sine = z0 + amp * np.sin(k * x)

for i in range(nb, nx+nb):
  for j in range(nb, nz+nb):
    if i >= sine[j-nb]:
      model[i, j] = 2000.0
    else:
      model[i, j] = 1400.0

# topo e base no mesmo loop
# preencher a laterais no mesmo loop

for j in range(nb, nx + nb):
    for i in range(nb):
        # top
        model[i, j] = model[nb, j]
        # bottom
        model[nz + nb + i, j] = model[nz + nb - 1, j]

for i in range(nzz):
    for j in range(nb):
        # left
        model[i, j] = model[i, nb]
        # right
        model[i, nx + nb + j] = model[i, nx + nb - 1]

plt.imshow(model, aspect="auto")

plt.tight_layout()
plt.show()

