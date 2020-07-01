module SegwayDyn

using OrdinaryDiffEq

#= State vector should be
x[1] : x position coordinate
x[2] : y position coordinate
x[3] : angle around z axis (negative yaw angle)
x[4] : angle around y axis (pitch angle)
x[5] : angular velocity around y axis
x[6] : forward velocity
x[7] : angular velocity around z axis
Input vector:
u[1] : normalised (in [-1, 1]) input voltaege to left motor
u[2] : normalised (in [-1, 1]) input voltaege to right motor
=#

include("body.jl")
include("drive.jl")
include("segway.jl")
export Segway
export segway_rk4, segway_rk4!

include("disc_sim.jl")
export dxdt_segway, segway_odestep

end # module SegwayRoboter
