# Functions pertaining to the curvilinear (spherical) RIC frame and the
#corresponding HCW equations

"""
    sph = ric2sph(tgtR, tgtRdot, chaseR, ric)

Transform from a linear RIC relative state to curvilinear spherical

Uses Schaub's formulation to go from RIC to Spherical

# Arguments
- `tgtR::AbstractFloat`: The norm of the target spacecraft position
- `tgtRdot::AbstractFloat`: The radial component of the target spacecraft velocity
- `chaseR::AbstractFloat`: The norm of the chase spacecraft position
- `ric::AbstractArray`: A 6-vector describing the RIC position and velocity of the chase
vehicle relative to the target
"""
function ric2sph(
    tgtR::AbstractFloat,
    tgtRdot::AbstractFloat,
    chaseR::AbstractFloat,
    ric::AbstractArray
)

    sph = zeros(6)

    sph[1] = chaseR - tgtR
    sph[2] = atan(ric[2], tgtR + ric[1])
    sph[3] = asin(ric[3] / chaseR)

    sph[4] = ((tgtR + ric[1]) * (tgtRdot + ric[4]) + ric[2] * ric[5] + ric[3] * ric[6])
    sph[4] /= chaseR
    sph[4] -= tgtRdot

    sph[5] = ((tgtR + ric[1]) * ric[5] - ric[2] * (tgtRdot + ric[4]))
    sph[5] /= ((tgtR + ric[1])^2 + ric[2]^2)

    sph[6] = ((tgtR + sph[1]) * ric[6] - ric[3] * (tgtRdot + sph[4]))
    sph[6] /= ((chaseR)^2 * sqrt(1 - (ric[3]^2 / ((chaseR)^2))))

    return sph
end

"""
    relativestate = eci2sph(target, chase)

Transform from ECI states to the chase's relative state about target in spherical
curvilinear coordinates

Uses Schaub's formation of the transformation, first converting to RIC then into spherical

# Arguments
- `target::AbstractArray`: A 6-vector describing the ECI position and velocity of the
target vehicle
- `chase::AbstractArray`: A 6-vector describing the ECI position and velocity of the
chase vehicle
"""
function eci2sph(target::AbstractArray, chase::AbstractArray)
    tgtR = norm(target[1:3])
    tgtRdot = dot(target[4:6], (target[1:3] / tgtR))
    chaseR = norm(chase[1:3])
    ric = eci2ric(target, chase)

    return ric2sph(tgtR, tgtRdot, chaseR, ric)
end


"""
    ric = sph2ric(tgtR, chaseR, sph)

Transform from a curvilinear spherical relative state to linear RIC

Uses Schaub's formulation to go from Spherical to RIC

# Arguments
- `tgtR::AbstractFloat`: The norm of the target spacecraft position
- `tgtRdot::AbstractFloat`: The radial component of the target spacecraft velocity
- `sph::AbstractArray`: A 6-vector describing the curvilinear spherical position and
velocity of the chase vehicle relative to the target
"""
function sph2ric(tgtR::AbstractFloat, tgtRdot::AbstractFloat, sph::AbstractArray)
    ric = zeros(6)
    chaseR = tgtR + sph[1]
    chaseRdot = tgtRdot + sph[4]

    ric[1] = chaseR * cos(sph[2]) * cos(sph[3]) - tgtR
    ric[2] = chaseR * sin(sph[2]) * cos(sph[3])
    ric[3] = chaseR * sin(sph[3])

    ric[4] = chaseRdot * cos(sph[2]) * cos(sph[3])
    ric[4] -= sph[5] * chaseR * sin(sph[2]) * cos(sph[3])
    ric[4] -= sph[6] * chaseR * cos(sph[2]) * sin(sph[3]) + tgtRdot

    ric[5] = chaseRdot * sin(sph[2]) * cos(sph[3])
    ric[5] += sph[5] * chaseR * cos(sph[2]) * cos(sph[3])
    ric[5] -= sph[6] * chaseR * sin(sph[2]) * sin(sph[3])

    ric[6] = chaseRdot * sin(sph[3]) + sph[6] * chaseR * cos(sph[3])

    return ric
end

"""
    chase = sph2eci(target, sph)

Transform from the chase's relative state about target in spherical curvilinear coordinates
to the chase's ECI state

Uses Schaub's formation of the transformation, first converting to RIC then into inertial

# Arguments
- `target::AbstractArray`: A 6-vector describing the ECI position and velocity of the
target vehicle
- `sph::AbstractArray`: A 6-vector describing the spherical relative position and velocity
of the chase vehicle
"""
function sph2eci(target::AbstractArray, sph::AbstractArray)
    tgtR = norm(target[1:3])
    tgtRdot = dot(target[4:6], (target[1:3] / tgtR))

    ric = sph2ric(tgtR, tgtRdot, sph)

    return ric2eci(target, ric)
end

"""
    phi = sphstm(dt, sma = 42164.15405)

Compute and return the state transition matrix for the Clohessy Wiltshire equations for use
on a curvilinear state with angular In-track and cross-track components.
"""
function sphstm(dt::AbstractFloat; sma::AbstractFloat=42164.15405)
    w = sqrt(3.986e5 / sma^3)
    ct = cos(w * dt)
    st = sin(w * dt)

    Phi = zeros(6, 6)
    Phi[1, 1] = 4 - 3 * ct
    Phi[1, 4] = st / w
    Phi[1, 5] = (1 - ct) * 2 * sma / w

    Phi[2, 1] = 6 * (st - w * dt) / sma
    Phi[2, 2] = 1
    Phi[2, 4] = (ct - 1) * 2 / w / sma
    Phi[2, 5] = 4 * st / w - 3 * dt

    Phi[3, 3] = ct
    Phi[3, 6] = st / w

    Phi[4, 1] = 3 * w * st
    Phi[4, 4] = ct
    Phi[4, 5] = 2 * st * sma

    Phi[5, 1] = 6 * w * (ct - 1) / sma
    Phi[5, 4] = -2 * st / sma
    Phi[5, 5] = 4 * ct - 3

    Phi[6, 3] = -w * st
    Phi[6, 6] = ct

    return Phi
end

"""
    dx = sphHill(x, p, t)

Compute and return the relative state derivatives for the ODE version of HCW

# Arguments
- `x::AbstractArray`: The relative position and velocity of the chase vehicle
- `p`: Tuple containing the target's semi-major axis in [1] and any applied accelerations in the RIC frame in [2]
- `t`: Time, unused within the function but required for ODEProblem
"""
function sphHill(x, p, t)
    mu = 398600.4418
    sma = p[1]
    u = p[2]

    n = sqrt(mu / sma^3)
    A = [
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
        3*n^2 0.0 0.0 0.0 2*n*sma 0.0
        0.0 0.0 0.0 -2*n/sma 0.0 0.0
        0.0 0.0 -n^2 0.0 0.0 0.0
    ]

    B = [
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0/sma 0.0
        0.0 0.0 0.0 0.0 0.0 1.0/sma
    ]'

    return A * x + B * u
end

"""
    newstate = HCW_Sph_Prop(relativestate, dt, sma, accel=[0, 0, 0])

Propagate a spherical relativestate forward in time by dt utilizing the HCW equations

# Arguments
- `relativestate::AbstractArray`: The relative spherical position and velocity of the
chase vehicle
- `dt::AbstractFloat`: The amount of time, in seconds, to propagate forward
- `sma::AbstractFloat`: The semi-major axis, in kilometers, of the target spacecraft
- `R0dot::AbstractFloat`: The target's velocity projected onto its position vector
- `accel::AbstractArray = [0, 0, 0]`: Any thrusting acceleration applied to the chase
vehicle over the timestep dt
- `flag::Bool = true`: If true uses the closed form version, otherwise ODE version
"""
function HCW_Sph_Prop(
    relativestate::AbstractArray,
    dt::AbstractFloat,
    sma::AbstractFloat;
    accel::AbstractArray=[0, 0, 0],
    flag::Bool=true,
)
    
    if flag
        Phi = sphstm(dt, sma=sma)
        return Phi * relativestate + dt * [[0, 0, 0]; accel]
    else
        t = (0, dt)
        prob = ODEProblem(sphHill, relativestate, t, (sma, accel))
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        return sol.u[end]
    end  
end


"""
    HCW_Sph_Prop!(relativestate, dt, sma, accel=[0, 0, 0])

Propagate a spherical relativestate forward in time by dt utilizing the HCW equations

# Arguments
- `relativestate::AbstractArray`: The relative spherical position and velocity of the
chase vehicle
- `dt::AbstractFloat`: The amount of time, in seconds, to propagate forward
- `sma::AbstractFloat`: The semi-major axis, in kilometers, of the target spacecraft
- `R0dot::AbstractFloat`: The target's velocity projected onto its position vector
- `accel::AbstractArray = [0, 0, 0]`: Any thrusting acceleration applied to the chase
vehicle over the timestep dt
- `flag::Bool = true`: If true uses the closed form version, otherwise ODE version
"""
function HCW_Sph_Prop!(
    relativestate::AbstractArray,
    dt::AbstractFloat,
    sma::AbstractFloat;
    accel::AbstractArray=[0, 0, 0]
)
    
    if flag
        phi = sphstm(dt, sma=sma)
        relativestate = phi * relativestate + dt * [[0, 0, 0]; accel]
    else
        t = (0, dt)
        prob = ODEProblem(sphHill, relativestate, t, (sma, accel))
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        relativestate = copy(sol.u[end])
    end  
end

function HCW_Sph_Solve(
    r0::AbstractArray,
    r1::AbstractArray,
    dt::AbstractFloat;
    sma::AbstractFloat=42164.15405
)

    phi = sphstm(dt, sma=sma)

    phiRR = phi[1:3, 1:3]
    phiRV = phi[1:3, 4:6]

    return phiRV \ (r1 - phiRR * r0)
end

"""
    meas = sphAzElMeasurements(state, tgtR, tgtRdot)

Generate azimuth and elevation measurements for a given RIC state

# Arguments
- `state::AbstractArray`: The relative position and velocity of the chase vehicle
- `tgtR::AbstractFloat`: The norm of the target spacecraft position
- `tgtRdot::AbstractFloat`: The radial component of the target spacecraft velocity
"""
function sphAzElMeasurements(
    state::AbstractVector, 
    tgtR::AbstractFloat, 
    tgtRdot::AbstractFloat,
)
    ric = sph2ric(tgtR, tgtRdot, state)

    return azElMeasurements(ric)
end
