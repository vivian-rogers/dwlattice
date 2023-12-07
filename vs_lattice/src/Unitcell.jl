# Specific imports.
include("DW_Lattice.jl")

# units
μm = 10^-6; GHz = 10^9; nm = 10^-9; #Ms_CoFeB = 1440.0*10^3; #A/m
Ms_CoFeB = 1273E3 # A/m
#=
#testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
#H_AFM_racetrack = constructHamiltonian(testSystem,2)

#AFM2R_lattice = init(; n_racetracks::Int, racetrack_positions::Vector{Float64}, orientations::Vector{Int},
#        a::Float64, ω₀::Vector{Float64}, PBC::Bool, w_RT::Float64, Ms::Float64, t::Float64, \alpha  w_DW::Float64=40*nm)
=#

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0,100]*nm, orientations=[1,-1], 
    a=200*nm, ω₀=0.0*[1.0,5.0]*GHz, PBC=true, w_RT=25*nm, Ms=Ms_CoFeB, t=3*nm, α=0.0015, w_DW=25.0*nm);


FM2R_lattice = init(n_racetracks=2, racetrack_positions=[0,50]*nm, orientations=[1,1], 
    a=100*nm, ω₀=2π*[10.0,10.0]*GHz, PBC=true, w_RT=30*nm, Ms=Ms_CoFeB, t=3*nm, α=0.015, w_DW=5*nm);

FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*nm, orientations=[1], 
    a=100*nm, ω₀=2π*[10.00]*GHz, PBC=true, w_RT=30*nm, Ms=Ms_CoFeB, t=3*nm, α=0.015, w_DW=5*nm);


    #=
AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0,2.5]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.1*[1,1]*GHz, ω₀=0.0*[0,0]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm, );

AFM3R_lattice = init(n_racetracks=3, racetrack_positions=[0, 2.8, 5.0]*μm, orientations=[1,-1,1], 
    a=7.5*μm, γ=0.5*[1,0.5,1]*GHz, ω₀=10*[0,2.0,2.1]*GHz, PBC=true, C=400.0*GHz, R₀=2.5*μm);

FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*μm, orientations=[1], 
    a=5*μm, γ=0.2*[1]*GHz, ω₀=10*[2.0]*GHz, PBC=true, C=80.0*GHz, R₀=2.5*μm);
=#