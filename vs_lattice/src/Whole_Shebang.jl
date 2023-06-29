# we want to write a code to, when given a unit cell with a certain set of parameters, 
# perform all of these interesting analyses!
include("DW_Position.jl")
include("Plot.jl")
include("Unitcell.jl")

function wholeAnalysis(system::DWLattice, Ncells::Int, race_i::Int, race_j::Int, NNs::Int, x::Function, f₀::Float64, name::String)
    # Make directory for the unit lattice if doesn't exist.
    if(!isdir("vs_lattice/img/" * name))
        mkpath("vs_lattice/img/" * name)
    end

    # Magic Numbers
    Δt = 1/f₀
    ω = 2*π/f₀    # angular freq.
    step = 1
    Npts = 100
    ϵ = 0.001

    # Position y(t)
    y = DWPosition(system, Ncells, race_i, race_j, x, f₀)
    # Transfer Function T
    T = transferFunc(system, Ncells, race_i, race_j, 0.0, 0.0, step, NNs=NNs, eval=false)

    # get the band structure and DOS of the periodic structure
    getBands(system, NNs, true, title=name)

    # now get the Bode Plot regarding the transfer function of one site on another
    BodePlot(T, ϵ*GHz, 25*GHz, Npts, title=name)

    # now convolve the impulse response with the input DW pos x(t) to get y(t)
    Plot_TwoFs(x, y, Δt, title=name)
end

f₀ = 1*GHz
wholeAnalysis(AFM2R_lattice, 10, 10, 7, 200, genFiring_neuron_DWpos(f₀), f₀, "AFM2R_Triangle")
wholeAnalysis(FM1R_lattice, 12, 8, 7, 200, genFiring_neuron_DWpos(f₀), f₀, "FM1R_Triangle")