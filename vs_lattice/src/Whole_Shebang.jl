# we want to write a code to, when given a unit cell with a certain set of parameters, 
# perform all of these interesting analyses!

function wholeAnalysis(system::DWLattice, Ncells::Int, race_i::Int, race_j::Int, name::String)
    path = "analyse_"*name*"/"
    # get the band structure and DOS of the periodic structure

    # now get the Bode Plot regarding the transfer function of one site on another

    # now convolve the impulse response with the input DW pos x(t) to get y(t)
    
end