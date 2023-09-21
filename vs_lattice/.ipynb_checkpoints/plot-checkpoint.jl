include("./src/DW_Lattice.jl")

"""
Set testing/ plotting env.
"""

testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
H_AFM_racetrack = constructHamiltonian(testSystem,2)

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0, 2.6]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.5*[1,1]*GHz, ω₀=1.0*[9,8]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm);
AFM3R_lattice = init(n_racetracks=3, racetrack_positions=[0, 2.8, 5.0]*μm, orientations=[1,-1,1], 
    a=7.5*μm, γ=0.5*[1,0.5,1]*GHz, ω₀=10*[0,2.0,2.1]*GHz, PBC=true, C=400.0*GHz, R₀=2.5*μm);
FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*μm, orientations=[1], 
    a=5*μm, γ=0.5*[1]*GHz, ω₀=10*[2.0]*GHz, PBC=true, C=80.0*GHz, R₀=2.5*μm);

"""
Testing/ plotting
"""
# getBands(AFM2R_lattice,200,true)
# getBands(FM1R_lattice,2000,true)
# getBands(AFM3R_lattice,50,true)

# test with AFM2R_lattice resonant 
# print(realSpaceArray(AFM2R_lattice, 3))

T, ω_val, mag, ϕ = transferFunc(AFM2R_lattice, 5, 1, 10, 0.0000001*GHz, 15.0*GHz, 200)

# Original Plot
plot(ω_val, mag)
title!("Magnitude")
xlabel!("ω (GHz)")
ylabel!("Pure Magnitude")
savefig("img/Orginal.png")

# 2-axis plot | Magnitude & Phase
plot(ω_val, mag, color=:green, label="Magnitude", xlabel="ω (GHz)", ylabel="Magnitude", legend=:topleft)
plot!(twinx(), ω_val, ϕ, color=:red, label=:"Phase", ylabel="Phase (rad)")
savefig("img/Multi_axis.png")

# Log(magnitude) plot
plot(ω_val, mag)
plot!(yscale=:log10)
title!("log(magnitude)")
xlabel!("log(ω)")
ylabel!("dB")
savefig("img/Log(magnitude).png")