import matplotlib as plt
import numpy as np

# DWTSW BAND STRUCTURE
def bandPlot (fft2D, label):
    from matplotlib.colors import LogNorm

    plt.figure(figsize=(10,6))
    plt.imshow(fft2D, aspect='auto', origin='lower', cmap="inferno", norm=LogNorm())

    plt.ylabel("$f$ (Hz)")
    plt.xlabel("$k$ (1/m)")

    plt.savefig("plots/" + label + "_bandStruct.png")

def posVtime (pos, T, dt, label, small=False):
    n_racetrack = range(pos.shape[0])

    # s, m version
    # T_vals = np.linspace(0, T, int(T/dt)+1)[:pos.shape[1]]
    # resize_track = pos * 0.25e-9

    # ns, nm version
    T_vals = np.linspace(0, T, int(T/dt)+1)[:pos.shape[1]] * 1e9
    resize_track = pos * 0.25

    if (small == True):
        plt.figure(figsize=(5, 3), dpi=300)
    else:
        plt.figure(figsize=(10, 6), dpi=300)

    for nth in n_racetrack:
        plt.plot(T_vals, resize_track[nth], label=str(nth+1)+'th')

    # ns, nm version
    plt.xlabel("Time (ns)")
    plt.ylabel("DW Position (nm)")
    plt.savefig("plots/" + label + "_thin_pos_v_t.png")