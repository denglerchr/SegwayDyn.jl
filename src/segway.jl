struct Segway{T}
    body::Body{T} # Body with its parameters
    driveleft::Drive{T} # left drive (motor+friction)
    driveright::Drive{T} # right drive (motor+friction)
end


Segway(T::DataType = Float64) = Segway{T}(Body(T), Drive(T), Drive(T))
Segway{T}() where {T} = Segway(T)


"""
function dxdt_segway(x::AbstractVector, u::AbstractVector, seg::Segway{T} = Segway(T))
Input u[1] is normalized voltage in [-1, 1] of the LEFT motor
u[2] is normalized voltage for the RIGHT motor
"""
function dxdt_segway(x::AbstractVector, u::AbstractVector, seg::Segway{T} = Segway(Float64)) where {T<:Number}
    dxdt = similar(x)
    return dxdt_segway!(dxdt, x, u, seg)
end


function dxdt_segway!(dxdt::AbstractVector, x::AbstractVector, u::AbstractVector, seg::Segway{T} = Segway(Float64)) where {T<:Number}
    ulim = clamp.(u, -one(eltype(u)), one(eltype(u)))

    # calculate the wheel angular velocities, relative to the body
    v_left = x[6] - x[7] * seg.body.b
    v_right = x[6] + x[7] * seg.body.b
    omega_right = v_right / seg.body.R - x[5]
    omega_left = v_left / seg.body.R - x[5]

    # Calculate torque from the drives. u[1] is input to the left wheel
    i_l = ulim[1] - seg.driveleft.kw * omega_left
    i_r = ulim[2] - seg.driveright.kw * omega_right
    taul = seg.driveleft(i_l, omega_left)
    taur = seg.driveright(i_r, omega_right)
    tau = [taul, taur]

    # Get state derivatives
    dxdt_body!( dxdt, x, tau, seg.body)

    return dxdt
end
