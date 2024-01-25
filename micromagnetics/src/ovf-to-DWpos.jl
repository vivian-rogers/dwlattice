using LinearAlgebra
using Glob
using DataFrames
using CSV

function bisect_array(Mzvals::Vector{Float64}, i_lo::Int, i_hi::Int)
        i₀ = (i_lo + i_hi) / 2.0 
        if mod(i₀,1.0) ≈ 0.5 && i_hi - i_lo ≈ 1.0
                return (i_lo, i_hi)
        end
        i₀ = Int(round(i₀))
        # ascending case
        if (Mzvals[i_lo] < Mzvals[i_hi])
                if Mzvals[i₀] > 0
                        (i_lo, i_hi) = bisect_array(Mzvals, i_lo, i₀)
                else
                        (i_lo, i_hi) = bisect_array(Mzvals, i₀, i_hi)
                end
        else # descending case
                if Mzvals[i₀] < 0
                        (i_lo, i_hi) = bisect_array(Mzvals, i_lo, i₀)
                else
                        (i_lo, i_hi) = bisect_array(Mzvals, i₀, i_hi)
                end
        end
        return (i_lo, i_hi)
end

function find_zero(Mzvals::Vector{Float64})
        # first we will solve this in terms of indices
        nMz = size(Mzvals)[1]
        i_lo, i_hi = bisect_array(Mzvals, 1, nMz)
        M_lo = Mzvals[i_lo]; M_hi = Mzvals[i_hi]
        i₀ = i_hi - (i_hi - i_lo)/(M_hi-M_lo)*M_hi
        return x₀ = i₀/nMz # is the domain wall position!
end

function get_positions(file)
        prev_row_zeros = true
        DWpositions = Float64[]
        DW_sub_positions = Float64[]
        iy = 1
        for row in CSV.Rows(file)
                #process it and turn it into Floats
                Mzvals = [parse(Float64, string(value)) for value in row]
                #Mvals = Float64.(row)
                if(sum(Mzvals.^2) != 0) #this row has magnetic data on it
                        if prev_row_zeros
                                iy = 1
                                DW_sub_positions = Float64[]
                                append!(DW_sub_positions, find_zero(Mzvals))
                        else
                                iy += 1
                                append!(DW_sub_positions, find_zero(Mzvals))
                        end
                        prev_row_zeros = false
                else
                        if !prev_row_zeros # we had data in the last row
                                append!(DWpositions, sum(DW_sub_positions)/iy)
                                #println("Got DW # $(size(DWpositions)[1])")
                        end
                        prev_row_zeros = true
                end
        end
        if !prev_row_zeros
                append!(DWpositions, sum(DW_sub_positions)/iy)
        end
        return DWpositions
end

function ovfs_to_DW_pos_csv(L::Float64=40.0, Δt::Float64=250E-10, nracetracks::Int=64)
	#run(`mumax3-convert -comp 2 -csv *.ovf`)
	csv_files = glob("./micromagnetics/testing/*.csv")
        csv_files = sort(csv_files)
        nT = size(csv_files)[1]
        # Initialize alldata as an array of arrays

        alldata = Matrix{Float64}(undef, nT, nracetracks+1)
        println("Processing with $(Threads.nthreads()) threads!")
        @Threads.threads for ifile in eachindex(csv_files)
            file = csv_files[ifile]
            println("Processing $file")
            DW_positions = L*get_positions(file)
            newrow = vcat(Δt*[ifile], DW_positions)
            # Ensure newrow is treated as a separate row with individual columns
            alldata[ifile,:] = newrow
        end
        
        # Convert alldata (array of arrays) into a matrix and then transpose it
        # Convert the transposed matrix into a DataFrame
        df = DataFrame(alldata, :auto)
        
        println("Printing DW positions and timestamps to DW_positions.csv!")
        # Save the DataFrame as a CSV file without writing headers
        CSV.write("./micromagnetics/testing/DW_positions.csv", df, writeheader=false)
        println("Done!")
end

ovfs_to_DW_pos_csv()
