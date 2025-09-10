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
nx, nz = 1150, 648
nxx, nzz = nx + 2 * nb, nz + 2 * nb

pSnapshots, vxSnapshots, vzSnapshots = getSnapshotNames(PATH)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))

ani = FuncAnimation(fig, update, frames=len(vxSnapshots), blit=False, interval=100)
ani.save("animation.mp4", writer="ffmpeg", fps=20)
plt.show()
