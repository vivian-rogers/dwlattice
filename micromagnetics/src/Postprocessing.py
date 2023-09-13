# GENERAL IMPORTS.
import numpy as np

# SPECIFIC IMPORTS.
import mumax as mx
import positionFinder as pf
import figures as fig

# PATH PARAMETERS
dataPath = "/work/09424/jiwoooop/ls6/mumax3-tacc-install/outputs/alp0temp150.mx3-general-mumax3.sl-13-39_Tuesday_August_29_2023/alp0temp150.out"

# NAME OF PROCESSED DATA
label = "a(0)_temp(150)"

# NUMERICAL PARAMETERS
fmax = 2e9       # maximum frequency (in Hz) of the sinc pulse
T    = 400e-9      # simulation time (longer -> better frequency resolution)
dt   = 100e-12     # the sample time

# Extract fields
outputdir  = dataPath
mx.ovf_to_npy(outputdir)
fields = mx.npy_to_dict(outputdir)

# POSITION OF DW IN LIST
# DW POSITION PARSING
pos = []
posE = []
posO = []

for key in sorted(fields.keys()):
    if "regions" in str(key):
        continue
    pos.append(pf.findDW(fields[key]))
    posE.append(pf.findDW_A(fields[key]))
    print("Processing: " + str(key))
    posO.append(pf.findDW_B(fields[key]))

# Saving copy
f = open("posData/pos_" + label + ".txt", 'w')
f.write(str(pos))

f = open("posData/posE_" + label + ".txt", 'w')
f.write(str(posE))

f = open("posData/posO_" + label + ".txt", 'w')
f.write(str(posO))


# Pos data extraction
from ast import literal_eval

# 1. IMPORT .TXT AS STRING
posStr = open("posData/pos_" + label + ".txt", "r").read()
posEvenStr = open("posData/posE_" + label + ".txt", 'r').read()
posOddStr = open("posData/posO_" + label + ".txt", "r").read()

# 2. CONVERT .TXT TO LIST
posList = literal_eval(posStr)
posEList = literal_eval(posEvenStr)
posOList = literal_eval(posOddStr)

# 3. CONVERT LIST TO NP ARRAY
pos = np.array(posList)
posE = np.array(posEList)
posO = np.array(posOList)

# Apply the two dimensional FFT
# DWTSW Band Structure, saves figure to plots dir
posE_fft = np.fft.fft2(posE)
posO_fft = np.fft.fft2(posO)
A = np.add(np.abs(posE_fft)**2, np.abs(posO_fft)**2)
fig.bandPlot(A, label)

# Plot positions along racetracks
by_track = np.transpose(pos)
fig.posVtime(by_track, T, dt, label)
fig.posVtime(by_track, T, dt, label, False)