module SemsFDTD

export SimData
export approx_stepsize
export new_sim
export timestep

"""
    SimData

    Struct containing the matrices required to run the FDTD simulation
"""
struct SimData
    e_x::AbstractArray{Float64,3}
    e_y::AbstractArray{Float64,3}
    e_z::AbstractArray{Float64,3}
    # H field definitions
    h_x::AbstractArray{Float64,3}
    h_y::AbstractArray{Float64,3}
    h_z::AbstractArray{Float64,3}
    # E field coefficients
    coeff_e_xe::AbstractArray{Float64,3}
    coeff_e_xh::AbstractArray{Float64,3}
    coeff_e_ye::AbstractArray{Float64,3}
    coeff_e_yh::AbstractArray{Float64,3}
    coeff_e_ze::AbstractArray{Float64,3}
    coeff_e_zh::AbstractArray{Float64,3}
    # H field coefficients
    coeff_h_xe::AbstractArray{Float64,3}
    coeff_h_xh::AbstractArray{Float64,3}
    coeff_h_ye::AbstractArray{Float64,3}
    coeff_h_yh::AbstractArray{Float64,3}
    coeff_h_ze::AbstractArray{Float64,3}
    coeff_h_zh::AbstractArray{Float64,3}
    # ABC X axis Planes
    abc_xp0_ty::AbstractArray{Float64,2}
    abc_xp0_tz::AbstractArray{Float64,2}
    abc_xp1_ty::AbstractArray{Float64,2}
    abc_xp1_tz::AbstractArray{Float64,2}
    # ABC Y axis Planes
    abc_yp0_tx::AbstractArray{Float64,2}
    abc_yp0_tz::AbstractArray{Float64,2}
    abc_yp1_tx::AbstractArray{Float64,2}
    abc_yp1_tz::AbstractArray{Float64,2}
    # ABC Z axis planes
    abc_zp0_tx::AbstractArray{Float64,2}
    abc_zp0_ty::AbstractArray{Float64,2}
    abc_zp1_tx::AbstractArray{Float64,2}
    abc_zp1_ty::AbstractArray{Float64,2}
end

# freespace impedance
vacuumZ = 376.730314 # ohms

# Courant limit in 3d space is 1/sqrt(3)
courant_limit = 1.0 / sqrt( 3.0 );

# first order ABC coefficient
abc_coeff = ( courant_limit - 1.0 ) / ( courant_limit + 1.0 )

# speed of light
c = 299792458

"""
    approx_stepsize - Give the approximate discrete spatial and time step sizes
    for a given minimum wavelength
"""
function approx_stepsize( min_wavelength::Float64 )
    delta_s = min_wavelength / 20
    delta_t = ( delta_s * courant_limit ) / c;
    return ( delta_s, delta_t )
end

"""
    new_sim - construct and initialize a new simulation assuming
"""
function new_sim( dim_x::Int, dim_y::Int, dim_z::Int )::SimData
    newsim = SimData( zeros( dim_x - 1, dim_y, dim_z ), # e_x
                      zeros( dim_x, dim_y - 1, dim_z ), # e_y
                      zeros( dim_x, dim_y, dim_z - 1 ), # e_z
                      zeros( dim_x, dim_y - 1, dim_z - 1 ), # h_x
                      zeros( dim_x - 1, dim_y, dim_z - 1 ), # h_y
                      zeros( dim_x - 1, dim_y - 1, dim_z ), # h_z
                      zeros( dim_x - 1, dim_y, dim_z ), # coeff_e_xe
                      zeros( dim_x - 1, dim_y, dim_z ), # coeff_e_xh
                      zeros( dim_x, dim_y - 1, dim_z ), # coeff_e_ye
                      zeros( dim_x, dim_y - 1, dim_z ), # coeff_e_yh
                      zeros( dim_x, dim_y, dim_z - 1 ), # coeff_e_ze
                      zeros( dim_x, dim_y, dim_z - 1 ), # coeff_e_zh
                      zeros( dim_x, dim_y - 1, dim_z - 1 ), # coeff_h_xe
                      zeros( dim_x, dim_y - 1, dim_z - 1 ), # coeff_h_xh
                      zeros( dim_x - 1, dim_y, dim_z - 1 ), # coeff_h_ye
                      zeros( dim_x - 1, dim_y, dim_z - 1 ), # coeff_h_yh
                      zeros( dim_x - 1, dim_y - 1, dim_z ), # coeff_h_ze
                      zeros( dim_x - 1, dim_y - 1, dim_z ), # coeff_h_zh
                      zeros( dim_y - 1, dim_z ), # abc_xp0_ty
                      zeros( dim_y, dim_z - 1 ), # abc_xp0_tz
                      zeros( dim_y - 1, dim_z ), # abc_xp1_ty
                      zeros( dim_y, dim_z - 1 ), # abc_xp1_tz
                      zeros( dim_x - 1, dim_z ), # abc_yp0_tx
                      zeros( dim_x, dim_z - 1 ), # abc_yp0_tz
                      zeros( dim_x - 1, dim_z ), # abc_yp1_tx
                      zeros( dim_x, dim_z - 1 ), # abc_yp1_tz
                      zeros( dim_x - 1, dim_y ), # abc_zp0_tx
                      zeros( dim_x, dim_y - 1 ), # abc_zp0_ty
                      zeros( dim_x - 1, dim_y ), # abc_zp1_tx
                      zeros( dim_x, dim_y - 1 ) ) # abc_zp1_ty

    # initialize e-field coefficients to free space for now
    newsim.coeff_e_xe .= 1.0
    newsim.coeff_e_xh .= courant_limit * vacuumZ

    return newsim
end

"""
    update_e_field -
"""
function update_e_field( sim::SimData )
    # x component
    sim.e_x[ :, begin + 1 : end - 1, begin + 1 : end - 1 ] =
        ( sim.coeff_e_xe[ :, begin + 1 : end - 1, begin + 1 : end - 1 ] .*
          sim.e_x[ :, begin + 1 : end - 1, begin + 1 : end - 1 ] ) .+
        sim.coeff_e_xh[ :, begin + 1 : end - 1, begin + 1 : end - 1 ] .*
        ( ( sim.h_z[ :, begin + 1 : end, begin + 1 : end - 1 ]
         .- sim.h_z[ :, begin : end - 1, begin + 1 : end - 1 ] )
        .- ( sim.h_y[ :, begin + 1 : end - 1, begin + 1 : end ]
          .- sim.h_y[ :, begin + 1 : end - 1, begin : end - 1 ] ) )
    # y component
    sim.e_y[ begin + 1 : end - 1, :, begin + 1 : end - 1 ] =
        ( sim.coeff_e_ye[ begin + 1 : end - 1, :, begin + 1 : end - 1 ] .*
          sim.e_y[ begin + 1 : end - 1, :, begin + 1 : end - 1 ] ) .+
        sim.coeff_e_yh[ begin + 1 : end - 1, :, begin + 1 : end - 1 ] .*
        ( ( sim.h_x[ begin + 1 : end - 1, :, begin + 1 : end ]
         .- sim.h_x[ begin + 1 : end - 1, :, begin : end - 1 ] )
        .- ( sim.h_z[ begin + 1 : end, :, begin + 1 : end - 1 ]
          .- sim.h_z[ begin : end - 1, :, begin + 1 : end - 1 ] ) )
    # z component
    sim.e_z[ begin + 1 : end - 1, begin + 1 : end - 1, : ] =
        ( sim.coeff_e_ze[ begin + 1 : end - 1, begin + 1 : end - 1, : ] .*
          sim.e_z[ begin + 1 : end - 1, begin + 1 : end - 1, : ] ) .+
        sim.coeff_e_zh[  begin + 1 : end - 1, begin + 1 : end - 1, : ] .*
        ( ( sim.h_y[ begin + 1 : end, begin + 1 : end - 1, : ]
         .- sim.h_y[ begin : end - 1, begin + 1 : end - 1, : ] )
        .- ( sim.h_x[ begin + 1 : end - 1, begin + 1 : end, : ]
          .- sim.h_x[ begin + 1 : end - 1, begin : end - 1, : ] ) )

    return nothing
end

"""
    update_h_field -
"""
function update_h_field( sim::SimData )
    # x component
    sim.h_x[:] = ( sim.coeff_h_xh .* sim.h_x ) .+ sim.coeff_h_xe .*
        ( ( sim.e_y[ :, :, begin + 1 : end ] .- sim.e_y[ :, :, begin : end - 1 ] )
        .- ( sim.e_z[ :, begin + 1 : end, : ] .- sim.e_z[ :, begin : end - 1, : ] ) )
    # y component
    sim.h_y[:] = ( sim.coeff_h_yh .* sim.h_y ) .+ sim.coeff_h_ye .*
        ( ( sim.e_z[ begin + 1 : end, :, : ] .- sim.e_z[ begin : end - 1, :, : ] )
        .- ( sim.e_x[ :, :, begin + 1 : end ] .- sim.e_x[ :, :, begin : end - 1 ] ) )
    # z component
    sim.h_y[:] = ( sim.coeff_h_zh .* sim.h_z ) .+ sim.coeff_h_ze .*
        ( ( sim.e_x[ :, begin + 1 : end, : ] .- sim.e_x[ :, begin : end - 1, : ] )
        .- ( sim.e_y[ begin + 1 : end, :, : ] .- sim.e_y[ begin : end - 1, :, : ] ) )

    return nothing
end

"""
    enforce_abc - enforce first order ABC at edges of box
"""
function enforce_abc( sim::SimData )
    # x plane 0
    sim.e_y[ begin, :, : ] = sim.abc_xp0_ty .+ abc_coeff * ( sim.e_y[ begin + 1, :, : ] - sim.e_y[ begin, :, : ] )
    sim.abc_xp0_ty[:] = sim.e_y[ begin + 1, :, : ]
    sim.e_z[ begin, :, : ] = sim.abc_xp0_tz .+ abc_coeff * ( sim.e_z[ begin + 1, :, : ] - sim.e_z[ begin, :, : ] )
    sim.abc_xp0_tz[:] = sim.e_z[ begin + 1, :, : ]
    # x plane 1
    sim.e_y[ end, :, : ] = sim.abc_xp1_ty .+ abc_coeff * ( sim.e_y[ end - 1, :, : ] - sim.e_y[ end, :, : ] )
    sim.abc_xp1_ty[:] = sim.e_y[ end - 1, :, : ]
    sim.e_z[ end, :, : ] = sim.abc_xp1_tz .+ abc_coeff * ( sim.e_z[ end - 1, :, : ] - sim.e_z[ end, :, : ] )
    sim.abc_xp1_tz[:] = sim.e_z[ end - 1, :, : ]
    # y plane 0
    sim.e_x[ :, begin, : ] = sim.abc_yp0_tx .+ abc_coeff * ( sim.e_x[ :, begin + 1, : ] - sim.e_x[ :, begin, : ] )
    sim.abc_yp0_tx[:] = sim.e_x[ :, begin + 1, : ]
    sim.e_z[ :, begin, : ] = sim.abc_yp0_tz .+ abc_coeff * ( sim.e_z[ :, begin + 1, : ] - sim.e_z[ :, begin, : ] )
    sim.abc_yp0_tz[:] = sim.e_z[ :, begin + 1, : ]
    # y plane 1
    sim.e_x[ :, end, : ] = sim.abc_yp1_tx .+ abc_coeff * ( sim.e_x[ :, end - 1, : ] - sim.e_x[ :, end, : ] )
    sim.abc_yp1_tx[:] = sim.e_x[ :, end - 1, : ]
    sim.e_z[ :, end, : ] = sim.abc_yp1_tz .+ abc_coeff * ( sim.e_z[ :, end - 1, : ] - sim.e_z[ :, end, : ] )
    sim.abc_yp1_tz[:] = sim.e_z[ :, end - 1, : ]
    # z plane 0
    sim.e_x[ :, :, begin ] = sim.abc_zp0_tx .+ abc_coeff * ( sim.e_x[ :, :, begin + 1 ] - sim.e_x[ :, :, begin ] )
    sim.abc_zp0_tx[:] = sim.e_x[ :, :, begin + 1 ]
    sim.e_y[ :, :, begin ] = sim.abc_zp0_ty .+ abc_coeff * ( sim.e_y[ :, :, begin + 1 ] - sim.e_y[ :, :, begin ] )
    sim.abc_zp0_ty[:] = sim.e_y[ :, :, begin + 1 ]
    # z plane 1
    sim.e_x[ :, :, end ] = sim.abc_zp1_tx .+ abc_coeff * ( sim.e_x[ :, :, end - 1 ] - sim.e_x[ :, :, end ] )
    sim.abc_zp1_tx[:] = sim.e_x[ :, :, end - 1 ]
    sim.e_y[ :, :, end ] = sim.abc_zp1_ty .+ abc_coeff * ( sim.e_y[ :, :, end - 1 ] - sim.e_y[ :, :, end ] )
    sim.abc_zp1_ty[:] = sim.e_y[ :, :, end - 1 ]

    return nothing
end

"""
    timestep - advance the simulation one timestep
"""
function timestep( sim::SimData )
    # enforce_pec()

    update_h_field( sim )
    update_e_field( sim )

    # handle any field dumping here

    # handle any sources here - TODO: before/after abc?

    enforce_abc( sim )

    return nothing
end

"""
    enforce_pec - apply zero electric field (PEC) to specified regions
        region argument expected in form of [ [ 1:2, 2:3, 3:4 ], [ 5:6, 7:8, 9:10 ] ]
"""
function enforce_pec( sim::SimData, regions::Vector{Vector{UnitRange{Int64}}} )
    for region in regions
        sim.e_x( region ) .= 0.0
        sim.e_y( region ) .= 0.0
        sim.e_z( region ) .= 0.0
    end

    return nothing
end

"""
    apply_source_node -
        TODO: does this make sense as a function
"""
function apply_source_node()
    return nothing
end

"""
    apply_array_nodes -
        TODO: does this make sense as a function
"""
function apply_array_nodes()
    return nothing
end

# function __init__()
#     print( "test" )
# 
#     tst = new_sim( 400, 400, 400 )
# end

end # module
