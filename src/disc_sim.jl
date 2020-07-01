# This file contains everything related to simulating in discrete time

function correct_fallen!(xt1, xt0)
    if -pi/2>=xt1[4]
        xt1[1:3] .= xt0[1:3]
        xt1[4] = -pi/2
        xt1[5:7] .= 0.0
    elseif pi/2<=xt1[4]
        xt1[1:3] .= xt0[1:3]
        xt1[4] = pi/2
        xt1[5:7] .= 0.0
    end
    return xt1
end

function correct_fallen2!(x, dx)
    if -pi/2>=x[4]+dx[4]
        x[4] = -pi/2
        x[5:7] .= 0.0
    elseif pi/2<=x[4]+dx[4]
        x[4] = pi/2
        x[5:7] .= 0.0
    else
        x .+= dx
    end
    return x
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
    return correct_fallen2!(xt, temp)
end


"""
function segway_rk4(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T

Do a discrete Runge Kutta 4 step, with step size dt. Returns next state without overwriting xt.
"""
function segway_rk4(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64)) where T
    xt2 = copy(xt)
    return segway_rk4!(xt2, ut, dt, seg)
end


# Stuff used for OrdinaryDifferentialEquations
function jacobian_fd(fun, x::AbstractVector{T}) where T
    jac = Array{T}(undef, length(x), length(x))
    return jacobian_fd!(jac, fun, x)
end

function jacobian_fd!(jac::AbstractArray{T}, fun, x::AbstractVector{T}) where T
    dx = similar(x)
    for i = 1:length(x)
        delta = max(T(1e-16), 10*eps(x[i]))
        fill!(dx, zero(T))
        dx[i] = delta
        jac[:, i] .= (fun(x+dx).-fun(x-dx))./(delta*2)
    end
    return jac
end


function segway_odestep(xt::AbstractVector, ut::AbstractVector, dt::Number, seg::Segway{T} = Segway(Float64); alg = OrdinaryDiffEq.VCABM()) where T
    # Adapt to notation of OrdinaryDifferentialEquations
    f = (dx, x, p, t)->dxdt_segway!(dx, x, ut, seg)
    jac = (J, x, p, t)->jacobian_fd!(J, x2->dxdt_segway(x2, ut, seg), x)
    jp = Array{T}(undef, length(xt), length(xt))

    fun = ODEFunction(f; jac=jac, jac_prototype = jp)
    prob = OrdinaryDiffEq.ODEProblem{true}(fun, xt, (zero(dt), dt))

    sol = OrdinaryDiffEq.solve(prob, alg; dt = dt)
    xt1 = sol.u[end]
    return correct_fallen!(xt1, xt)
end
