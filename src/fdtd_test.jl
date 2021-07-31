
include( "SemsFDTD.jl" )

testsim = SemsFDTD.new_sim( 400, 400, 400 )

SemsFDTD.timestep( testsim )

@time begin
    for i = 1:25
        SemsFDTD.timestep( testsim )
        print( ' ', i )
    end
end
