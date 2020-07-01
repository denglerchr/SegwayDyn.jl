"""
Drive struct, contains parameters for motor constant and friction.
Friction is modeled as visous friction. Motor constant is assumed not constant, to account for static friction.
(There is a bend in the characteristic curve).
"""
struct Drive{S}
    km::S # motor torque constant
    kw::S # Counter emf
    cfric1::S #friction aprameter
    cfric2::S #friction aprameter
    cfric3::S #friction aprameter
end

Drive(S::DataType; km = 0.61, kw = 0.018131147540983605, cfric1 = 0.24, cfric2 = 2.0, cfric3 = 0.4) = Drive{S}(S(km), S(kw), S(cfric1), S(cfric2), S(cfric3))
Drive{S}(; args...) where {S}= Drive(S; args...)


"""
Input is normalised voltage [-1, 1] and rotation speed of the wheel, output is torque on the body
"""
function (drive::Drive{S})(i::Number, omega::Number) where {S}
    # Constant * current
    tau = drive.km*i

    # Friction (Coulomb+viscous)
    tau -= drive.cfric1*tanh(drive.cfric2*omega)*exp(-drive.cfric3*omega^2)
    return tau
end
