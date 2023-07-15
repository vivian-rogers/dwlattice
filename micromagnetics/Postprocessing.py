# imports
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import Rbf

'''
Functions
'''
def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""

    from pandas import read_table
    
    table = read_table(filename)
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table


def read_mumax3_ovffiles(outputdir):
    """Load all ovffiles in outputdir into a dictionary of numpy arrays 
    with the ovffilename (without extension) as key"""
    
    from subprocess import run, PIPE, STDOUT
    from glob import glob
    from os import path
    from numpy import load

    # convert all ovf files in the output directory to numpy files
    p = run(["mumax3-convert","-numpy",outputdir+"/*.ovf"], stdout=PIPE, stderr=STDOUT)
    if p.returncode != 0:
        print(p.stdout.decode('UTF-8'))

    # read the numpy files (the converted ovf files)
    fields = {}
    for npyfile in glob(outputdir+"/*.npy"):
        key = path.splitext(path.basename(npyfile))[0]
        fields[key] = load(npyfile)
    
    return fields


'''
Extract table and fields
'''
name = "dwlattice/micromagnetics/DW_Lattice"

scriptfile = name + ".txt" 
outputdir  = name + ".out"

table = read_mumax3_table(outputdir + "/table.txt") 
fields = read_mumax3_ovffiles(outputdir)


'''
Return the value of DW position in a cell.
(By finding a point where mz == 0)
'''
def bisect(f, x_lo, x_hi, n):
    # midpoint
    x_o = (x_lo + x_hi) / 2.0
    # base case, run enough to converge.
    if n == 0:
        return x_o
    # base case, at x_o is 0
    if f(x_o) == 0:
        return x_o
    
    # ascending (-1 --> 1)
    if f(x_lo) < f(x_hi):
        if f(x_o) > 0:
            x_o = bisect(f, x_lo, x_o, n-1)
        else:
            x_o = bisect(f, x_o, x_hi, n-1)
    # descending (1 --> -1)
    else:
        if f(x_o) < 0:
            x_o = bisect(f, x_lo, x_o, n-1)
        else:
            x_o = bisect(f, x_o, x_hi, n-1)
    return x_o

'''
Return array of DW positions within each racetrack.
@params     m (field), i.e. fields["m000000"]
'''
def findDW(m):
    # *Global* values/ variables
    x_total = m.shape[3]
    y_total = m.shape[2]
    x_mid = int(m.shape[3]/2)

    DW_array = []
    DW_pos = []

    nDW = 0.0
    flag = False

    for iy in range(y_total):
        mx = m[0,0,iy,x_mid]
        my = m[1,0,iy,x_mid]
        mz = m[2,0,iy,x_mid]
        empty = (mx==0) and (my==0) and (mz==0)

        if empty:
            # Just got off the racetrack.
            if flag == True:
                DW_array.append(sum(DW_pos) / nDW)
                DW_pos = []
                nDW = 0
                flag = False
        else:
            # Just entered the racetrack.
            if flag == False:
                flag = True
            f = Rbf(np.arange(x_total), m[2,0,iy,:]) 
            DW_pos.append(bisect(f, 0, x_total, 15))
            nDW += 1
    
    return DW_array

def findDW_A(m):
    # *Global* values/ variables
    x_total = m.shape[3]
    y_total = m.shape[2]
    x_mid = int(m.shape[3]/2)

    DW_array = []
    DW_pos = []

    nDW = 0.0
    flag = False
    count = -1

    for iy in range(y_total):
        mx = m[0,0,iy,x_mid]
        my = m[1,0,iy,x_mid]
        mz = m[2,0,iy,x_mid]
        empty = (mx==0) and (my==0) and (mz==0)

        if empty:
            # Just got off the racetrack.
            if flag == True:
                if count%2 == 0:
                    DW_array.append(sum(DW_pos) / nDW)
                DW_pos = []
                nDW = 0
                flag = False
        else:
            # Just entered the racetrack.
            if flag == False:
                flag = True
                count += 1
            f = Rbf(np.arange(x_total), m[2,0,iy,:]) 
            DW_pos.append(bisect(f, 0, x_total, 15))
            nDW += 1
    
    return DW_array

def findDW_B(m):
    # *Global* values/ variables
    x_total = m.shape[3]
    y_total = m.shape[2]
    x_mid = int(m.shape[3]/2)

    DW_array = []
    DW_pos = []

    nDW = 0.0
    flag = False
    count = -1

    for iy in range(y_total):
        mx = m[0,0,iy,x_mid]
        my = m[1,0,iy,x_mid]
        mz = m[2,0,iy,x_mid]
        empty = (mx==0) and (my==0) and (mz==0)

        if empty:
            # Just got off the racetrack.
            if flag == True:
                if count%2 != 0:
                    DW_array.append(sum(DW_pos) / nDW)
                DW_pos = []
                nDW = 0
                flag = False
        else:
            # Just entered the racetrack.
            if flag == False:
                flag = True
                count += 1
            f = Rbf(np.arange(x_total), m[2,0,iy,:]) 
            DW_pos.append(bisect(f, 0, x_total, 15))
            nDW += 1
    
    return DW_array
'''
DW Position Parsing
'''
# Stack all snapshots of the magnetization on top of each other
# snaps = len(fields.keys())
# posA = np.stack([findDW_A(fields[key]) for key in iter.islice(sorted(fields.keys()), snaps-1)])
# posB = np.stack([findDW_B(fields[key]) for key in iter.islice(sorted(fields.keys()), snaps-1)])

pos = []
posA = []
posB = []

for key in sorted(fields.keys()):
    if key != "regions000000":
        pos.append(findDW(fields[key]))
        posA.append(findDW_A(fields[key]))
        print("Processing: " + str(key))
        posB.append(findDW_B(fields[key]))

posA_fft = np.fft.fft2(posA)
posB_fft = np.fft.fft2(posB)
posA_fft = np.fft.fftshift(posA_fft)
posB_fft = np.fft.fftshift(posB_fft)

A = np.add(np.abs(posA_fft)**2, np.abs(posB_fft)**2)
# A = np.abs(posA_fft)**2
# B = np.abs(posB_fft)**2


'''
DWTSW Band Structure
'''
from matplotlib.colors import LogNorm

plt.figure(figsize=(10,6))

# Show the intensity plot of the 2D FFT
# extent = [ -(2*np.pi)/(2*dx), (2*np.pi)/(2*dx), -1/(2*dt), 1/(2*dt)] # extent of k values and frequencies

# plt.imshow(np.abs(posA_fft)**2, extent=extent, aspect='auto', origin='lower', cmap="inferno")
# plt.imshow(np.abs(posB_fft)**2, extent=extent, aspect='auto', origin='lower', cmap="inferno")
plt.imshow(A, aspect='auto', origin='lower', cmap="inferno", norm=LogNorm())

# plt.xlim([-2e8,2e8])
# plt.ylim([0,fmax])
plt.ylabel("$f$ (Hz)")
plt.xlabel("$k$ (1/m)")

plt.show()

# Plot positions along racetracks
by_track = np.transpose(pos)

'''
Position v. Time
'''
n_racetrack = range(by_track.shape[0])      # 16
T_vals = np.linspace(0, T, int(T/dt)+1)

plt.figure(figsize=(10, 6), dpi=300)

for nth in n_racetrack:
    plt.plot(T_vals, by_track[nth], label=str(nth+1)+'th')

plt.legend(loc="upper left", fontsize='5')