#Can use this to linearisze the model, in case one wants to design a linear controller
using SegwayDyn
using ForwardDiff

seg = Segway(Float64)
A = ForwardDiff.jacobian(x->dxdt_segway(x, zeros(2), seg), zeros(7))
@show A

#B = ForwardDiff.jacobian(u->dxdt_segway(zeros(7), u, seg), zeros(2)) #some type problem

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

B = myjacobian(u->dxdt_segway(zeros(7), u, seg), zeros(2))
@show B
nothing
