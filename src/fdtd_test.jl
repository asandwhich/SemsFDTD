
using Plots
using WriteVTK
gr()

include( "SemsFDTD.jl" )

testsim = SemsFDTD.new_sim( Float32, 1000, 1000, 25 )

SemsFDTD.timestep( testsim )

@time begin
    for i = 1:50
        SemsFDTD.timestep( testsim )
        print( ' ', i )
    end
end

heatmap( testsim.coeff_e_xe[:,:,2], aspect_ratio = 1 )
# vtk_write_array( "test_file_out", rand( 100, 100, 100 )
