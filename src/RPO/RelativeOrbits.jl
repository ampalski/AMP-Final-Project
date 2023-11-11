module RelativeOrbits

using LinearAlgebra: norm, cross, dot
using DifferentialEquations

μ = 3.986e5
rate(a) = sqrt(μ / a^3)

# Wrap an unbounded angle to the interval [0, 2pi]
function _wrapTo2pi(input::Float64)
    p2 = 2 * pi
    while input < 0
        input += p2
    end

    return input % p2
end

export stm, sphstm
export Hill, sphHill
export HCW_Prop, HCW_Prop!
export HCW_Sph_Prop, HCW_Sph_Prop!
export HCW_Solve, HCW_Sph_Solve
export eci2ric, ric2eci
export ric2sph, sph2ric
export eci2sph, sph2eci

include("Rectilinear.jl")
include("Spherical.jl")

end # module RelativeOrbits