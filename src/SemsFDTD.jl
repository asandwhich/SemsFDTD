module SemsFDTD

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


function new_sim( dim_x::UInt32, dim_y::UInt32, dim_z::UInt32 )::SimData
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

    # free space impedance
    fsi = 
    # initialize e-field coefficients to free space for now
    newsim.coeff_e_xe .= 1.0
    newsim.coeff_e_xh = courantlimit * fsi
end

"""
    update_e_field -
"""
function update_e_field()
end

"""
    update_h_field -
"""
function update_h_field()
end

"""
    enforce_abc -
"""
function enforce_abc()
end

"""
    enforce_pec -
"""
function enforce_pec()
end

"""
    apply_source_node -
"""
function apply_source_node()
end

"""
    apply_array_nodes -
"""
function apply_array_nodes()
end

function __init__()
    print( "test" )

    tst = new_sim( 400, 400, 400 )
end

end # module
