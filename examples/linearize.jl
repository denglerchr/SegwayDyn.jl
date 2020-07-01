#Can use this to linearisze the model, in case one wants to design a linear controller
using SegwayDyn

seg = SegwayDyn.Segway(Float64)

function myjacobian(fun, x)
    deriv = Array{Vector}(undef, 0)
    dx = similar(x)
    for i = 1:length(x)
        delta = max(1e-16, 10*eps(x[i]))
        fill!(dx, 0)
        dx[i] = delta
        push!(deriv, (fun(x+dx).-fun(x-dx))./(2*delta))
    end
    return hcat(deriv...)
end

A = myjacobian(x->SegwayDyn.dxdt_segway(x, zeros(2), seg), zeros(9))
B = myjacobian(u->SegwayDyn.dxdt_segway(zeros(9), u, seg), zeros(2))
nothing
