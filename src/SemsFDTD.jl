module SemsFDTD


"""
    NetworkData

    Struct containing information describing a microwave network
"""
struct NetworkData
    data::AbstractArray{ComplexF64,3}
    freq::AbstractArray{Float64,1}
    paramType::NetParameterType
    ref_imp::Float64
    noiseData::AbstractArray{Float64,2}
end


"""
"""
struct SimData
    e_x::AbstractArray{Float64,3}
    e_y::AbstractArray{Float64,3}
    e_z;
    # H field definitions
    h_x;
    h_y;
    h_z;
    # E field coefficients
    coeff_e_xe;
    coeff_e_xh;
    coeff_e_ye;
    coeff_e_yh;
    coeff_e_ze;
    coeff_e_zh;
    # H field coefficients
    coeff_h_xe;
    coeff_h_xh;
    coeff_h_ye;
    coeff_h_yh;
    coeff_h_ze;
    coeff_h_zh;
    # ABC X axis Planes
    abc_xp0_ty;
    abc_xp0_tz;
    abc_xp1_ty;
    abc_xp1_tz;
    # ABC Y axis Planes
    abc_yp0_tx;
    abc_yp0_tz;
    abc_yp1_tx;
    abc_yp1_tz;
    # ABC Z axis planes
    abc_zp0_tx;
    abc_zp0_ty;
    abc_zp1_tx;
    abc_zp1_ty;
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
end

end # module
