# Functions pertaining to the linear RIC frame and the corresponding HCW equations


"""
    phi = stm(dt, sma = 42164.15405)

Compute and return the state transition matrix for the Clohessy Wiltshire equations.
"""
function stm(dt::AbstractFloat; sma::AbstractFloat=42164.15405)
    w = sqrt(3.986e5 / sma^3)
    ct = cos(w * dt)
    st = sin(w * dt)

    Phi = zeros(6, 6)
    Phi[1, 1] = 4 - 3 * ct
    Phi[1, 4] = st / w
    Phi[1, 5] = (1 - ct) * 2 / w

    Phi[2, 1] = 6 * (st - w * dt)
    Phi[2, 2] = 1
    Phi[2, 4] = (ct - 1) * 2 / w
    Phi[2, 5] = 4 * st / w - 3 * dt

    Phi[3, 3] = ct
    Phi[3, 6] = st / w

    Phi[4, 1] = 3 * w * st
    Phi[4, 4] = ct
    Phi[4, 5] = 2 * st

    Phi[5, 1] = 6 * w * (ct - 1)
    Phi[5, 4] = -2 * st
    Phi[5, 5] = 4 * ct - 3

    Phi[6, 3] = -w * st
    Phi[6, 6] = ct

    return Phi
end

"""
    dx = Hill(x, p, t)

Compute and return the relative state derivatives for the ODE version of HCW

# Arguments
- `x::AbstractArray`: The relative position and velocity of the chase vehicle
- `p`: Tuple containing the target's semi-major axis in [1] and any applied accelerations in the RIC frame in [2]
- `t`: Time, unused within the function but required for ODEProblem
"""
function Hill(x::AbstractVector, p, t)
    mu = 398600.4418
    sma = p[1]
    u = p[2]

    n = sqrt(mu / sma^3)
    A = [
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
        3*n^2 0.0 0.0 0.0 2*n 0.0
        0.0 0.0 0.0 -2*n 0.0 0.0
        0.0 0.0 -n^2 0.0 0.0 0.0
    ]

    B = [
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
    ]'

    return A * x + B * u
end

"""
    newstate = HCW_Prop(relativestate, dt, sma, accel=[0, 0, 0])

Propagate relativestate forward in time by dt utilizing the HCW equations

# Arguments
- `relativestate::AbstractArray`: The relative position and velocity of the chase vehicle
- `dt::AbstractFloat`: The amount of time, in seconds, to propagate forward
- `sma::AbstractFloat`: The semi-major axis, in kilometers, of the target spacecraft
- `accel::AbstractArray = [0, 0, 0]`: Any thrusting acceleration applied to the chase vehicle over the timestep dt
- `flag::Bool = true`: If true uses the closed form version, otherwise ODE version
"""
function HCW_Prop(
    relativestate::AbstractArray,
    dt::AbstractFloat,
    sma::AbstractFloat;
    accel::AbstractArray=[0, 0, 0],
    flag::Bool=true,
)

    if flag
        phi = stm(dt, sma=sma)
        return phi * relativestate + dt * [[0, 0, 0]; accel]
    else
        t = (0, dt)
        prob = ODEProblem(Recti_Hill, relativestate, t, (sma, accel))
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        return sol.u[end]
    end
end

"""
    HCW_Prop!(relativestate, dt, sma, accel = [0, 0, 0])

Propagate relativestate forward in time by dt utilizing the HCW equations

# Arguments
- `relativestate::AbstractArray`: The relative position and velocity of the chase vehicle
- `dt::AbstractFloat`: The amount of time, in seconds, to propagate forward
- `sma::AbstractFloat`: The semi-major axis, in kilometers, of the target spacecraft
- `accel = [0, 0, 0]`: Any thrusting acceleration applied to the chase vehicle over the
timestep dt
- `flag::Bool = true`: If true uses the closed form version, otherwise ODE version
"""
function HCW_Prop!(
    relativestate::AbstractArray,
    dt::AbstractFloat,
    sma::AbstractFloat;
    accel=[0, 0, 0],
    flag=true,
)
    if flag
        phi = stm(dt, sma=sma)
        relativestate = phi * relativestate + dt * [[0, 0, 0]; accel]
    else
        t = (0, dt)
        prob = ODEProblem(Recti_Hill, relativestate, t, (sma, accel))
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        relativestate = copy(sol.u[end])
    end

end

"""
    relativestate = eci2ric(target, chase)

Transform from ECI states to the chase's relative state about target in RIC

# Arguments
- `target::AbstractArray`: A 6-vector describing the ECI position and velocity of the
target vehicle
- `chase::AbstractArray`: A 6-vector describing the ECI position and velocity of the
chase vehicle
"""
function eci2ric(target::AbstractArray, chase::AbstractArray)
    R = unit(target[1:3])
    C = unit(cross(R, target[4:6]))
    I = unit(cross(C, R))

    posdiff = chase[1:3] - target[1:3]
    veldiff = chase[4:6] - target[4:6]

    ricpos = [R'; I'; C'] * posdiff

    # Velocity corrects for rotational component through transport theorem

    h_vec = cross(target[1:3], target[4:6])
    h = norm(h_vec)
    fdot = h / norm(target[1:3])^2
    omega = [0, 0, fdot]

    relvel = [R'; I'; C'] * veldiff - cross(omega, ricpos)

    return [ricpos; relvel]
end

"""
    chase = ric2eci(target, relativestate)

Transform from the chase's relative state about target in RIC to its ECI state

# Arguments
- `target::AbstractArray`: A 6-vector describing the ECI position and velocity of the
target vehicle
- `relativestate::AbstractArray`: A 6-vector describing the RIC relative position and
velocity of the chase vehicle
"""
function ric2eci(target::AbstractArray, relativestate::AbstractArray)
    eci = copy(target)

    R = unit(target[1:3])
    C = unit(cross(target[1:3], target[4:6]))
    I = unit(cross(C, R))

    RIC = [R I C]

    posdiff = RIC * relativestate[1:3]

    eci[1:3] += posdiff

    h_vec = cross(target[1:3], target[4:6])
    h = norm(h_vec)
    fdot = h / norm(target[1:3])^2
    omega = [0, 0, fdot]

    relvel = RIC * (relativestate[4:6] + cross(omega, relativestate[1:3]))

    eci[4:6] += relvel
    return eci
end

"""
    v0 = HCW_Solve(r0, r1, dt, sma = 42164.15405)

Solves the velocity required to get from position r0 to position r1 in dt seconds.

# Arguments
- `r0::AbstractArray`: A 3-vector describing the starting RIC position
- `r1::AbstractArray`: A 3-vector describing the ending RIC position
- `dt::Float64`: The propagation time, in seconds
- `sma::Float64`: The semi-major axis of the target vehicle
"""
function HCW_Solve(
    r0::AbstractArray,
    r1::AbstractArray,
    dt::Float64;
    sma::Float64=42164.15405
)

    phi = stm(dt, sma=sma)

    phiRR = phi[1:3, 1:3]
    phiRV = phi[1:3, 4:6]

    v0 = phiRV \ (r1 - phiRR * r0)

    return v0

end

"""
    meas = azElMeasurements(state)

Generate azimuth and elevation measurements for a given RIC state

# Arguments
- `state::AbstractArray`: The relative position and velocity of the chase vehicle
"""
function azElMeasurements(state::AbstractVector)
    measurement = zeros(2)

    measurement[1] = atan(state[2], state[1])
    measurement[2] = atan(state[3], sqrt(state[1]^2 + state[2]^2))

    return _wrapTo2pi.(measurement)
end