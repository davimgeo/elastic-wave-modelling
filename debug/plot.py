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
nrec = 172;

seismogram = np.fromfile(
    "data/output/seismogram_txx_1150x648.bin", dtype=np.float32, count=nt*nrec
).reshape([nt, nrec], order='F')

perc = 99

scale_min = np.percentile(seismogram, 100-perc)
scale_max = np.percentile(seismogram, perc)

TRACE = 80

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 8), sharey=True)

ax[0].imshow(seismogram, cmap="Greys", aspect="auto", vmin=scale_min, vmax=scale_max)
ax[0].plot(TRACE * np.ones(nt), np.arange(nt), 'r--')

ax[1].plot(seismogram[:, TRACE], np.arange(nt))

plt.tight_layout()
plt.show()


