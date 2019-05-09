struct Segway{T}
    body::Body{T} # Body with its parameters
    driveleft::Drive{T} # left drive (motor+friction)
    driveright::Drive{T} # right drive (motor+friction)
end


Segway(T::DataType) = Segway{T}(Body(T), Drive(T), Drive(T))
Segway{T}() where {T} = Segway(T)


function dxdt_segway(x::AbstractVector, u::AbstractVector, seg::Segway{T} = Segway(T)) where {T<:Number}
    out = similar(x)

    # calculate the wheel angular velocities
    v_left = x[6] - x[7] * seg.body.b
    v_right = x[6] + x[7] * seg.body.b
    omega_right = v_right / seg.body.R - x[5]
    omega_left = v_left / seg.body.R - x[5]

    # Calculate torque from the drives. u[1] is input tot he left wheel
    taul = seg.driveleft(u[1], omega_left)
    taur = seg.driveright(u[2], omega_right)
    tau = [taul, taur]

    # Get state derivatives
    return dxdt_body(x, tau, seg.body)
end


"""
function segway_rk4!(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T

Do a discrete Runge Kutta 4 step, with step size dt. Overwrites (and returns) xt, which now holds the next state.
"""
function segway_rk4!(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T
    temp = similar(xt)

    #get partial results dy1-dy4
    dy1 = dxdt_segway(xt, ut, seg);
    temp .= dy1.*(dt/2) .+ xt;

    dy2 = dxdt_segway(temp, ut, seg);
    temp .= dy2.*(dt/2) .+ xt;

    dy3 = dxdt_segway(temp, ut, seg);
    temp .= dy3.*dt .+ xt;

    dy4 = dxdt_segway(temp, ut, seg);

    # update xt
    @. temp = dt*(dy1+2.0*(dy2+dy3)+dy4)/6.0

    #Check if segway has fallen
    if -pi/2>=xt[4]+temp[4]
        xt[4] = -pi/2
        xt[5:7] .= 0.0
    elseif pi/2<=xt[4]+temp[4]
        xt[4] = pi/2
        xt[5:7] .= 0.0
    else
        xt .+= temp
    end
    return xt
end


"""
function segway_rk4(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T

Do a discrete Runge Kutta 4 step, with step size dt. Returns next state without overwriting xt.
"""
function segway_rk4(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T
    xt2 = copy(xt)
    return segway_rk4!(xt2, ut, dt, seg)
end
