using SegwayDyn, ForwardDiff
import SegwayDyn: Body, Drive
customtype = Float32

## Create a custom segway
body = Body(customtype)
driveleft = Drive(customtype; cfric = 0.004)
driveright = Drive(customtype; cfric = 0.00446)
segway = Segway(body, driveleft, driveright)

## Create a time discrete PID controller for angle and position

mutable struct PID{T}
    P::T
    I::T
    D::T
    int_temp::T
    diff_temp::T
end

PID{T}(p, i, d) where {T} = PID(T(p), T(i), T(d), zero(T), zero(T))

function (pid::PID)(x::Number, dt::Number)
    pid.int_temp += dt*x
    diffx = (x-pid.diff_temp)/dt
    pid.diff_temp = x
    return pid.P * x + pid.I * pid.int_temp + pid.D * diffx
end

pidv = PID{customtype}(0.22, 0.0, 0.001)
pidalpha = PID{customtype}(2.0, 0.01, 0.0)

## Simulate
dt = customtype(0.01)
t = 0:dt:10
X = zeros(customtype, 7, length(t))
x0 = zeros(customtype, 7); x0[4] = 2*pi/180
X[:, 1] .= x0
U = zeros(customtype, 2, length(t))
for i = 2:length(t)
    U[:, i-1] .= clamp(pidalpha(X[4, i-1], dt) + pidv(X[6, i-1], dt), -one(customtype), one(customtype))
    X[:, i] .= segway_rk4(X[:, i-1], U[:, i-1], dt, segway)
end

## Plot results
using Plots
fig = plot(layout = (2, 1))
plot!(fig[1], X[[4, 6], :]')
plot!(fig[2], U')
