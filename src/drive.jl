"""
Drive struct, contains parameters for motor constant and friction.
Friction is modeled as visous friction. Motor constant is assumed not constant, to account for static friction.
(There is a bend in the characteristic curve).
"""
struct Drive{T}
    kl::T # motor constant at low input voltages
    kh::T # motor constant at higehr input voltages
    bend::T # "Voltage" at which motor "constant" changes (must be in [0, 1])
    cfric::T # Visouc friction + motor induction
end

Drive(T::DataType; kl = 0.138, kh = 0.465, bend = 0.06367, cfric = 0.0074) = Drive{T}(T(kl), T(kh), T(bend), T(cfric))
Drive{T}(; args...) where {T}= Drive(T; args...)


"""
Input is normalised voltage [-1, 1] and rotation speed of the wheel, output is torque on the body
"""
function (drive::Drive{T})(u::Number, omega::Number) where {T}
    out = zero(T)

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
