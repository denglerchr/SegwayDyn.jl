using SegwayDyn, ForwardDiff
customtype = Float32

## Create a custom segway
body = SegwayDyn.Body(customtype; bend = 0.2)
drivel = SegwayDyn.Drive(customtype; cfric = 0.004)
driver = SegwayDyn.Drive(customtype; cfric = 0.00446)
segway = Segway(body, drivel, driver)

## Create a stabilising controller
# Linearise model
A = ForwardDiff.J( x-> dxdt_segway( x , zeros(customtype, 2), segway), zeros(customtype, 7))

# Create LQR controller

## Simulate
