"""
Body struct, on case one wants to have custom parameters. Properties:
	Mb # kg # Masse des Körpers
	Mw # kg # Masse der Räder
	R # m #Radius der Räder
	cz # m #Schwerpunkt des Körpers Höhe
	b # m #Half the distance between the wheels (OO_w in paper)
	Ixx # kg m^2 #Inertia of the body
	Iyy # kg m^2 #Inertia of the body
	Izz # kg m^2 #Inertia of the body
	Iwa # kg m^2 #Inertia of the wheels about their axis
	Iwd # kg m^2 #Inertia of the wheels about a diameter
"""
struct Body{T}
    Mb::T # kg # Masse des Körpers
    Mw::T # kg # Masse der Räder
    R::T # m #Radius der Räder
    cz::T # m #Schwerpunkt des Körpers Höhe
    b::T # m #Half the distance between the wheels (OO_w in paper)
    Ixx::T # kg m^2 #Inertia of the body
    Iyy::T # kg m^2 #Inertia of the body
    Izz::T # kg m^2 #Inertia of the body
    Iwa::T # kg m^2 #Inertia of the wheels about their axis
    Iwd::T # kg m^2 #Inertia of the wheels about a diameter
end


Body{T}() where T = Body(T)


function Body(T::DataType = Float64;
	Mb = 1.76,
	Mw = 0.147,
	R = 0.07,
	cz = 0.077,
	b = 0.1985/2,
	Ixx = (0.166^2+0.21^2)*Mb/12 + Mb*(cz-0.07)^2,
	Iyy = (0.072^2+0.21^2)*Mb/12 + Mb*(cz-0.07)^2,
	Izz = (0.072^2+0.166^2)*Mb/12,
	Iwa = Mw*R^2/2,
    Iwd = Mw*b^2)

    return Body{T}(T(Mb), T(Mw), T(R), T(cz), T(b), T(Ixx), T(Iyy), T(Izz), T(Iwa), T(Iwd))
end

function Body(params::AbstractVector{T}) where{T}
    @assert(length(params) == 10)
    return Body{T}(params...)
end

"""
Body dynamics, taken from Pathak2005.
State vector is [x0, y0, phi, alpha, dalpha, v, dphi]
"""
@inline function dxdt_body(x::AbstractVector, u::AbstractVector, body::Body{T} = Body()) where {T<:Number}
    g = T(9.81) #N/kg

    # Rename variables
    phi = x[3]
    alpha = x[4]
    dalpha = x[5]
    v = x[6]
    dphi = x[7]

    # Rigid Body Dynamics
    Dalpha = body.Mb^2*cos(alpha)^2*body.cz^2*body.R^2+((-body.Mb^2-2*body.Mw*body.Mb)*body.cz^2-2*body.Iyy*body.Mw-body.Iyy*body.Mb)*body.R^2-2*body.Mb*body.cz^2*body.Iwa-2*body.Iyy*body.Iwa
    Galpha = (-body.Mb*body.cz^2+body.Izz-body.Ixx)*body.R^2*cos(alpha)^2 + (body.Mb*body.cz^2+body.Ixx+2*body.Iwd+2*body.b^2*body.Mw)*body.R^2+2*body.b^2*body.Iwa
    H = 1/2*body.Mb*body.R^2*body.Izz+body.Iwa*body.Izz-body.Mw*body.R^2*body.Ixx-body.Iwa*body.Ixx-body.Mb*body.cz^2*body.Mw*body.R^2 - body.Mb*body.cz^2*body.Iwa -1/2*body.Mb*body.R^2*body.Ixx+body.Mw*body.R^2*body.Izz
    Kalpha = (-4*body.Iyy*body.Mb*body.R^2*body.cz-3*body.R^2*body.Mb^2*body.cz^3+body.Mb*body.R^2*body.cz*(body.Ixx-body.Izz))*sin(alpha) + (body.Mb*body.R^2*body.cz*(body.Ixx-body.Izz)+body.R^2*body.Mb^2*body.cz^3)*sin(3*alpha)
    f21 = sin(2*alpha)*dphi^2*H/Dalpha + body.Mb^2*body.cz^2*body.R^2*sin(2*alpha)*dalpha^2/(2*Dalpha) + (-2*body.Mb^2*body.R^2*body.cz-4*body.Iwa*body.Mb*body.cz-4*body.Mw*body.R^2*body.Mb*body.cz)*g*sin(alpha)/(2*Dalpha)
    f22 = Kalpha*dphi^2+(body.Mb^2*body.cz^2*body.R^2*g*sin(2*alpha))/(2*Dalpha)+(-4*body.Iyy*body.Mb*body.R^2*body.cz-4*body.R^2*body.Mb^2*body.cz^3)*sin(alpha)*dalpha^2/(4*Dalpha)
    f23 = (-(body.Ixx-body.Izz)*body.R^2-body.Mb*body.cz^2*body.R^2)*sin(2*alpha)*dalpha*dphi/Galpha-sin(alpha)*body.R^2*body.Mb*body.cz*v*dphi/Galpha
    g21 = (u[1] + u[2])*(body.Mb*body.R^2+2*body.Mw*body.R^2+2*body.Iwa+body.Mb*cos(alpha)*body.cz*body.R)/Dalpha
    g22 = -(u[1] + u[2])*body.R*(body.Mb*cos(alpha)*body.cz*body.R+body.Iyy+body.Mb*body.cz^2)/Dalpha
    g23 = (u[1] - u[2])*body.R*body.b/Galpha

    xdot = similar(x)
    xdot[1] = cos(phi)*v
    xdot[2] = sin(phi)*v
    xdot[3] = x[7]
    xdot[4] = x[5]
    xdot[5] = f21 + g21
    xdot[6] = f22 + g22
    xdot[7] = f23 + g23

    return xdot
end

# Should not really be used
function body_rk4(xt::AbstractVector, ut::AbstractVector, dt::Number, body::Body = Body())
    xt1 = copy(xt)
    body_rk4!(xt1, ut, dt, body)
    return xt1
end


function body_rk4!(xt::AbstractVector, ut::AbstractVector, dt::Number, body::Body = Body())
    temp = similar(xt)

    #get partial results dy1-dy4
    dy1 = dxdt_body(xt, ut, body);
    temp .= dy1.*(dt/2) .+ xt;

    dy2 = dxdt_body(temp, ut, body);
    temp .= dy2.*(dt/2) .+ xt;

    dy3 = dxdt_body(temp, ut, body);
    temp .= dy3.*dt .+ xt;

    dy4 = dxdt_body(temp, ut, body);

    # update xt
    @. temp = dt*(dy1+2.0*(dy2+dy3)+dy4)/6.0

    #Check if segway has fallen #TODO limits anpassen
    if -pi/2>=xt[4]+temp[4]
        xt[4] = -pi/2
		xt[5:7] .= 0.0
    elseif pi/2<=xt[4]+temp[4]
        xt[4] = pi/2
		xt[5:7] .= 0.0
    else
        xt .+= temp
    end
    return nothing
end
