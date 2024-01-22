using LinearAlgebra
using Glob

function find_zero(Mzvals::Float64)
end

function get_positions(file)

function ovfs_to_DW_pos_csv()
	run(`mumax3-convert -csv *.ovf`)
	csv_files = glob("*.csv", "./")
	for file in csv_files
		println("Processing $file")
		DW_positions = get_positions(file)		
			
end
