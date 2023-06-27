# we want to write a code to, when given a unit cell with a certain set of parameters, 
# perform all of these interesting analyses!
include("DW_Position.jl")
include("Plot.jl")
include("Unitcell.jl")

function wholeAnalysis(system::DWLattice, Ncells::Int, race_i::Int, race_j::Int, x::Function, f₀::Float64, name::String)
    # Magic Numbers
    Δt = 1/f₀

    # Position y(t)
    y = DWPosition(system, Ncells, race_i, race_j, x, f₀)

    # get the band structure and DOS of the periodic structure

    # now get the Bode Plot regarding the transfer function of one site on another

    # now convolve the impulse response with the input DW pos x(t) to get y(t)
    Plot_TwoFs(x, y, Δt, title=name)
end

wholeAnalysis(FM1R_lattice, 2, 2, 1, genFiring_neuron_DWpos(f₀), 1.0, "Triangle")