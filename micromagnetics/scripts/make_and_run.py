import argparse
import os
import subprocess
import sys
import numpy as np

def main():
    
    # get the arguments in order:
    # need to 
    # module load intel
    # module load python3
    parser = argparse.ArgumentParser(description='Modify mumax3 input file parameters.')
    
    # janky supergarbage to set default material parameters
    material = sys.argv[1]
    parts = material.split('=')
    # The part after the '=' is the second element of the resulting list
    material = parts[1] if len(parts) > 1 else None

    # okay this is a bit silly. Let's do a first pass and get the material params if the material is defined
    #args = parser.parse_args()
    mat_params = {}
    print(material) 
    # do a jank match-case because python 3.9 doesn't support this
    if material == "YIG":
        # citations 10.1088/0022-3727/48/1/015001
        mat_params.update({"alpha":6.15E-5,  "Msat":141E3,       "Aex":3.7E-12,  "Ku1":0})
    else: # we just assume it's cofeb
        mat_params.update({"alpha":0.015,    "Msat": 1273E3,    "Aex": 10E-12,   "Ku1":1.5e6})
    print(mat_params) 
    parser.add_argument('--material', type=str, default="cofeb", help='Value for Exchange? (default: cofeb)')
    parser.add_argument('--alpha', type=float, default=mat_params["alpha"], help='Value for alpha (default: 0.01)')
    parser.add_argument('--Msat', type=float, default=mat_params["Msat"], help='Value for Saturation Magnetization (default: 1273E3)')
    parser.add_argument('--Aex', type=float, default=mat_params["Aex"], help='Value for Exchange stiffness J/m  (default: 10E-12)')
    parser.add_argument('--Ku1', type=float, default=mat_params["Ku1"], help='Value for Saturation Magnetization (default: 1.5E6)')
    parser.add_argument('--delta', type=float, default=30, help='Value for delta (nm) (default: 30)')
    parser.add_argument('--width', type=float, default=10, help='Value for width (nm) (default: 10)')
    parser.add_argument('--lattice_constant', type=float, default=60, help='Value for lattice_constant (nm) (default: 60)')
    parser.add_argument('--temp', type=int, default=273, help='Value for temp (default: 273)')

    # Parse the arguments
    args = parser.parse_args()

    for arg in vars(args):
        globals()[arg] = getattr(args, arg)


    Bz = 0.005*np.exp(-temp/100)
    dirname = f"{material}_delta#{delta}_width#{width}_latticeconst#{lattice_constant}_temp#{temp}" 

    fstring = f"""// NUMERICAL PARAMETERS
    fmax := 1.5e9    // maximum frequency (in Hz) of the sinc pulse
    T_tot := 400e-9  // simulation time (longer -> better frequency resolution)
    del_t := 150e-12 // the sample time

    // Racetrack crystal settings:
    length := 200e-9
    width := {width}e-9
    delta := {delta}e-9
    lattice_constant := {lattice_constant}e-9
    N_unitcells := 32
    Nx := 256
    Ny := 1024
    Nz := 1
    Lx := length
    Ly := N_unitcells * lattice_constant
    dX := length / Nx
    dY := Ly / Ny
    dZ := 3e-9 / Nz

    setGridSize(Nx, Ny, Nz)
    setCellSize(dX, dY, dZ)

    // Material Constants
    Bz := {Bz}
    Msat = {Msat}
    Aex = {Aex}
    anisU = vector(0, 0, 1)
    Ku1 = {Ku1}
    alpha = {alpha}
    Xi = 0.2
    demagAccuracy = 6
    pinning := false
    AFM := true

    two_racetrack := rect(length, delta+width)
    racetrack := rect(length, width)

    gap := delta - width

    singleR := racetrack.repeat(0, lattice_constant, 0).transl(0, width/2, 0)
    twoR := singleR.add(singleR.transl(0, delta, 0))

    if AFM {{
            setgeom(twoR)

            // Define regions
            defRegion(2, singleR.transl(0, delta, 0))
            save(regions)

            m.setRegion(2, twoDomain(0, 0, -1, 1, 1, 0, 0, 0, 1))
    }} else {{
            setgeom(singleR)
    }}

    defRegion(1, singleR)
    defRegion(3, racetrack.transl(0, -Ly/2+width/2, 0))
    save(regions)

    B_ext.setregion(3, vector(0, 0, Bz*sinc(2*pi*fmax*(t-T_tot/2))))

    m.setRegion(1, twoDomain(0, 0, 1, 1, 1, 0, 0, 0, -1))
    m.setRegion(3, twoDomain(0, 0, 1, 1, 1, 0, 0, 0, -1))

    setPBC(0, 2, 0)
    relax()

    // Schedule output & save results
    autosave(m, del_t)
    tableadd(e_total)
    tableautosave(del_t)

    // Run for 1ns with current through the sample
    temp.set({temp})
    j = vector(0, 0, 0)
    pol = 1
    run(T_tot)"""

    # Create a directory
    print(f"made dir {dirname}")
    os.makedirs(dirname, exist_ok=True)
    # Change the current working directory to the new directory
    os.chdir(dirname)
    with open('run.mx3', 'w') as file:
        file.write(fstring)
    
    command = "bash mumax3.sh racetrack-lattice.sl run.mx3"

    # Start the subprocess and get the output stream
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Print output in real-time
    while True:
        output = process.stdout.readline()
        if output == b'' and process.poll() is not None:
            break
        if output:
            print(output.strip().decode())

    process.poll()
    os.chdir("../")
    

main()
