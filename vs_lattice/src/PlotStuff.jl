function BodePlot(T::Function, startω::Float64, finalω::Float64, step::Int)
    ω_vals = LinRange(startω,finalω,step)
    mag_vals = (abs∘T).(ω_vals)
    ϕ_vals = (angle∘T).(ω_vals)
    #ω_val, mag_vals, ϕ_vals = transferFunc(system, n_lattice, race_i, race_j, startω, finalω, step); lw=4
    c_left=:black; c_right=:black; lw=4;
    plot(ω_vals, mag_vals, color=:blue, xlabel="ω (GHz)", label="Magnitude", legend_foreground_color=nothing,
        ylabel="|H(ω)| (dB)", lw=lw, legend=false, yscale=:log10, dpi=300, thickness_scaling = 1.5,
        y_foreground_color_text =c_left, y_foreground_color_border = c_left, y_foreground_color_axis= c_left, grid=false)
    plot!(twinx(), ω_vals, ϕ_vals, color=:red, alpha=0.6, label=:"Phase", legend=false, legend_foreground_color=nothing,
        ylabel="∠H(ω) (Rad)", lw=lw, yticks=([-π,0,π],["-π","0","π"]), grid=true, gridlinewidth=1, ylim=[-3.2,3.2],
        y_foreground_color_text =c_right, y_foreground_color_border = c_right, y_foreground_color_axis = c_right)
    savefig("./vs_lattice/img/BodePlot.png")
end

function BodePlot(system::DWLattice, n_lattice::Int, race_i::Int, race_j::Int, startω::Float64, finalω::Float64, step::Int, NNs::Int=50)
    ω_val, mag_vals, ϕ_vals = transferFunc(system, n_lattice, race_i, race_j, startω, finalω, step); lw=4
    c_left=:black; c_right=:black
    plot(ω_val, mag_vals, color=:blue, xlabel="ω (GHz)", label="Magnitude", legend_foreground_color=nothing,
        ylabel="|H(ω)| (dB)", lw=lw, legend=false, yscale=:log10, dpi=300, thickness_scaling = 1.5,
        y_foreground_color_text =c_left, y_foreground_color_border = c_left, y_foreground_color_axis= c_left, grid=false)
    plot!(twinx(), ω_val, ϕ_vals, color=:red, alpha=0.6, label=:"Phase", legend=false, legend_foreground_color=nothing,
        ylabel="∠H(ω) (Rad)", lw=lw, yticks=([-π,0,π],["-π","0","π"]), grid=true, gridlinewidth=1, ylim=[-3.2,3.2],
        y_foreground_color_text =c_right, y_foreground_color_border = c_right, y_foreground_color_axis = c_right)
end

function Plot_TwoFs(f1::Function, f2::Function,Δt::Float64,npts::Int=500)
    # Plot
    tvals = LinRange(0,Δt,npts)
    plot(tvals, f2, lw=4,c="",thickness_scaling=1.5)
    plot!(tvals, f1, lw=4,c="red", alpha=0.6)
    title!("x(t) and y(t) plotter")
    xlabel!("t (ns)")
    ylabel!("DW position (nm)")
    savefig("vs_lattice/img/plotio.png")
end


function Plot_Spectrum(harmonics, magnitude, title="")
    # Plot
    plot(harmonics, real.(magnitude), line=:stem,c="blue")
    plot!(harmonics, imag.(magnitude), line=:stem,c="red")
    title!("Magnitude of Harmonics")
    xlabel!("n (ᵗʰ harmonic)")
    ylabel!("Magnitude")
    savefig("vs_lattice/img/"*title*" FourierSeries.png")
end