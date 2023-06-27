# units
μm = 10^-6; GHz = 1; 

testSystem = DWLattice(2, [0, 2.5]*μm, [1,-1], 5*μm, [0.1,0.1]*GHz, [10,15]*GHz, true, 10^8, 2.5*μm);
H_AFM_racetrack = constructHamiltonian(testSystem,2)

AFM2R_lattice = init(n_racetracks=2, racetrack_positions=[0, 2.5]*μm, orientations=[1,-1], 
    a=5*μm, γ=0.1*[1,1]*GHz, ω₀=0.0*[9,8]*GHz, PBC=true, C=10.0*GHz, R₀=2.5*μm);
AFM3R_lattice = init(n_racetracks=3, racetrack_positions=[0, 2.8, 5.0]*μm, orientations=[1,-1,1], 
    a=7.5*μm, γ=0.5*[1,0.5,1]*GHz, ω₀=10*[0,2.0,2.1]*GHz, PBC=true, C=400.0*GHz, R₀=2.5*μm);
FM1R_lattice = init(n_racetracks=1, racetrack_positions=[0]*μm, orientations=[1], 
    a=5*μm, γ=0.5*[1]*GHz, ω₀=10*[2.0]*GHz, PBC=true, C=80.0*GHz, R₀=2.5*μm);