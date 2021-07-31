
include( "SemsFDTD.jl" )

testsim = SemsFDTD.new_sim( Float32, 400, 400, 400 )

SemsFDTD.timestep( testsim )

@time begin
    for i = 1:10
        SemsFDTD.timestep( testsim )
        print( ' ', i )
    end
end
