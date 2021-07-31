
include( "SemsFDTD.jl" )

testsim = SemsFDTD.new_sim( 60, 60, 60 )

SemsFDTD.timestep( testsim )

@time begin
    for i = 1:10
        SemsFDTD.timestep( testsim )
        print( ' ', i )
    end
end
