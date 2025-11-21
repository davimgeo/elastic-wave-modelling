import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def extract_tid(fname: str):
    return int(fname.split("tid_")[1].split(".")[0])

def getSnapshotNames(PATH: str):
    fields = {"p": [], "vx": [], "vz": []}
    for filename in os.listdir(PATH):
        matched = False
        for key in fields:
            if key in filename:
                fields[key].append(filename)
                matched = True
                break
        if not matched:
            raise ValueError(f"Invalid file inside snapshots: {filename}")
        
    for key in fields:
        fields[key].sort(key=extract_tid)
    return fields["p"], fields["vx"], fields["vz"]

def update(i):
    ax.clear()

    salt = np.fromfile(
        "data/output/vp.bin", dtype=np.float32, count=nxx*nzz
    ).reshape([nzz, nxx], order='F')
    
    snap = np.fromfile(
        os.path.join(PATH, vxSnapshots[i]), dtype=np.float32, count=nxx*nzz
    ).reshape([nzz, nxx], order='F')
    
    img_model = ax.imshow(salt, cmap='viridis', aspect='auto', alpha=0.5)
    img_snap = ax.imshow(snap, cmap='Greys', aspect='auto', alpha=0.7)

    ax.set_title(f"{pSnapshots[i]}")
    
    return img_model, img_snap,

PATH = "data/output/snapshots/"

nb = 100
nx, nz = 1701, 351
nxx, nzz = nx + 2 * nb, nz + 2 * nb

pSnapshots, vxSnapshots, vzSnapshots = getSnapshotNames(PATH)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))

ani = FuncAnimation(fig, update, frames=len(vxSnapshots), blit=False, interval=100)
#ani.save("animation.mp4", writer="ffmpeg", fps=20)
plt.show()

nt = 5001
nrec = 172
dt = 1e-3
offset = 10

seismogram = np.fromfile(
    "data/output/seismogram_txx_1150x648.bin",
    dtype=np.float32,
    count=nt * nrec
).reshape([nt, nrec], order='F')

perc = 99

scale_min = np.percentile(seismogram, 100 - perc)
scale_max = np.percentile(seismogram, perc)

fig, ax = plt.subplots(figsize=(10, 8))

tloc = np.linspace(0, nt - 1, 11, dtype = int)
tlab = np.around(tloc * dt, decimals = 1)

xloc = np.linspace(0, nrec - 1, 9)
xlab = np.array(offset * xloc, dtype = int)

scale_min = np.percentile(seismogram, 100 - perc)
scale_max = np.percentile(seismogram, perc)

img = ax.imshow(
    seismogram, aspect = "auto", cmap="Greys", 
    vmin=scale_min, vmax=scale_max
)

cbar = fig.colorbar(img, ax = ax, extend = 'neither')
cbar.minorticks_on()

ax.set_yticks(tloc)
ax.set_yticklabels(tlab)

ax.set_xticks(xloc)
ax.set_xticklabels(xlab)

ax.set_title("Seismogram", fontsize = 18)
ax.set_xlabel("Distance [m]", fontsize = 15)
ax.set_ylabel("Two Way Time [s]", fontsize = 15)

plt.tight_layout()
plt.show()



