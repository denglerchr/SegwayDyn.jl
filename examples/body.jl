using SegwayDyn

body = SegwayDyn.Body(Float32; cz = 0.72) # same as Body(Float32)

x0 = randn(Float32, 7)
u0 = randn(Float32, 2)
x1 = SegwayDyn.body_rk4(x0, u0, 0.01, body)
x1 = SegwayDyn.body_rk4(x0, u0, 0.01)
