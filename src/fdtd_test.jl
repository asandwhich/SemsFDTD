
using Plots
gr()

include( "SemsFDTD.jl" )

testsim = SemsFDTD.new_sim( Float16, 400, 400, 50 )

SemsFDTD.timestep( testsim )

@time begin
    for i = 1:5
        SemsFDTD.timestep( testsim )
        print( ' ', i )
    end
end

heatmap( testsim.coeff_e_xe[:,:,2], aspect_ratio = 1 )

