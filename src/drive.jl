"""
Drive struct, contains parameters for motor constant and friction.
Friction is modeled as visous friction. Motor constant is assumed not constant, to account for static friction.
(There is a bend in the characteristic curve).
"""
struct Drive{S}
    kl::S # motor constant at low input voltages
    kh::S # motor constant at higehr input voltages
    bend::S # "Voltage" at which motor "constant" changes (must be in [0, 1])
    cfric::S # Visouc friction + motor induction
    T::S # Timeconstant of PT1
end

Drive(S::DataType; kl = 0.138, kh = 0.465, bend = 0.06367, cfric = 0.0074, T = 0.01) = Drive{S}(S(kl), S(kh), S(bend), S(cfric), S(T))
Drive{S}(; args...) where {S}= Drive(S; args...)


"""
Input is normalised voltage [-1, 1] and rotation speed of the wheel, output is torque on the body
"""
function (drive::Drive{S})(u::Number, omega::Number) where {S}
    out = zero(S)

    # Motor constant
    if abs(u) < drive.bend
        out += drive.kl*u
    else
        out += sign(u)*( drive.kl*drive.bend + drive.kh*( abs(u) - drive.bend ) )
    end

    # Friction and motor induction
    out -= drive.cfric*omega
    return out
end
