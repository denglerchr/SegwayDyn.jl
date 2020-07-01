using SegwayDyn

testdrive = SegwayDyn.Drive(Float32; km = 0.2) # same as Body(Float32)

tau1 = testdrive(randn(), randn())
tau2 = testdrive(randn(Float32), randn(Float32))
