using StatsBase
include("DW_Lattice.jl")

"""
Bands
"""
#H_AFM_racetrack = constructHamiltonian(testSystem,2)

function plot1DBands(H::Function,system::DWLattice,nk::Int,Broadening::Bool=false, nE::Int = 200, DOS=true)
    a = system.a; nE = 400; nk = 400
    kvals = LinRange(-π/a,π/a,nk)
    ys = []
    xs = []
    maxE = maximum(imag.(eigvals(H(π/a)) ∪ eigvals(H(0.0))))
    if(Broadening)
        # get some reasonable maximum for the energies
        Evals = LinRange(0,maxE*1.1,nE)
        n = 2*system.n_racetracks[]
        bands = zeros(nk,nE)
        plottingBands = zeros(nk,nE)
        @Threads.threads for ik in eachindex(kvals)
            k = kvals[ik]
            Hatk = H(k)
            DOS_k(E) = (-1/π)imag.(tr(inv(E*I(n) - im*Hatk)))
            DOS_k_E = DOS_k.(Evals)
            bands[ik,:] = DOS_k_E
        end
        #bands = (1/maximum(bands))*bands
        #fig = heatmap(kvals*(a/π), Evals, bands', clims=(0,maximum(bands)), xlabel="k (π/a)", ylabel="Frequency (GHz)")
        fig = heatmap(kvals*(a/π), Evals, bands', clims=(0,quantile(vec(bands),0.995)), xlabel="k (π/a)", ylabel="ω (10⁹ rad/s)",dpi=300)
    end
    #@Threads.threads for k in kvals
    #    Es = imag(eigvals(H(k)))
    #    append!(ys,Es/GHz)
    #    append!(xs,k*ones(size(Es)))
    #end
    #fig = scatter!(xs*(a/π),ys, legend=false,ylims=(-0.001,maxE*1.1), xlims=(-1,1), c="white", markersize=2.0, markerstrokewidth=0)
    if(DOS)
        DOS_tot = [sum(bands[:,iE]) for iE = 1:nE]/nk; DOS_tot = DOS_tot
        fDOS = plot(DOS_tot,Evals,legend=false,xlabel="DOS (a.u.)", ylabel="ω (10⁹ rad/s)",ylims=(-0.0001,maxE*1.1),grid=false,xlims=(0,quantile(vec(DOS_tot),0.99)+0.1), lw=3,dpi=300)
    end
    fig2 = plot(fig,fDOS,layout=grid(1,2, widths=(5/8,3/8)), size=(800,300),margin=5mm,dpi=300)
    return fig2
end

function getBands(system::DWLattice, NNs::Int, broadening::Bool=true; title::String)
    H = constructHamiltonian(system,NNs)
    fig = plot1DBands(H,system,500,broadening,300)
    savefig("vs_lattice/img/"*title*"/_bands.png")
end

"""
Bode Plot
"""

function BodePlot(T::Function, startω::Float64, finalω::Float64, step::Int; title::String)
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
    savefig("./vs_lattice/img/"*title*"/_Bode.png")
end

function BodePlot(system::DWLattice, n_lattice::Int, race_i::Int, race_j::Int, startω::Float64, finalω::Float64, step::Int, NNs::Int=50; title::String)
    ω_val, mag_vals, ϕ_vals = transferFunc(system, n_lattice, race_i, race_j, startω, finalω, step); lw=4
    c_left=:black; c_right=:black
    plot(ω_val, mag_vals, color=:blue, xlabel="ω (GHz)", label="Magnitude", legend_foreground_color=nothing,
        ylabel="|H(ω)| (dB)", lw=lw, legend=false, yscale=:log10, dpi=300, thickness_scaling = 1.5,
        y_foreground_color_text =c_left, y_foreground_color_border = c_left, y_foreground_color_axis= c_left, grid=false)
    plot!(twinx(), ω_val, ϕ_vals, color=:red, alpha=0.6, label=:"Phase", legend=false, legend_foreground_color=nothing,
        ylabel="∠H(ω) (Rad)", lw=lw, yticks=([-π,0,π],["-π","0","π"]), grid=true, gridlinewidth=1, ylim=[-3.2,3.2],
        y_foreground_color_text =c_right, y_foreground_color_border = c_right, y_foreground_color_axis = c_right)
end

"""
Position y(t)
"""

function Plot_TwoFs(f1::Function, f2::Function, Δt::Float64, npts::Int=500; title::String)
    # Plot
    depin = 0.075
    tvals = LinRange(0,Δt,npts)
    f2vals = f2.(tvals)
    f1vals = f1.(tvals)
    f1vals .-= sum(f1vals)/size(f1vals)[1];
    plot(tvals, f1vals, lw=4,c="darkorchid", alpha=0.6, thickness_scaling=1.5,grid=false,dpi=300,legend=:topright,label="xᵢ(t)")
    plot!(tvals, f2vals, lw=2,c="green", legend_foreground_color=nothing,label="xⱼ(t)",
        xlim=(0.0,Δt),ylim=(-1.1,1.1), yticks=([-1,-depin,depin,1],["-L/2","-xₚ","xₚ","L/2"]))
    plot!([0,Δt],[depin,depin],lw=0.5,ls=:dash,c="black",label=nothing)
    plot!([0,Δt],[-depin,-depin],lw=0.5,ls=:dash,c="black",label=nothing)
    #title!("x(t) and y(t) plotter")
    xlabel!("t (ns)")
    ylabel!("DW position")
    savefig("vs_lattice/img/"*title*"/_y(t).png")
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
