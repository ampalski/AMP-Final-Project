#Need:
# DONE: struct that outlines what configuration the target is in, where the cones are, etc.
# DONE: instantaneous point check against those KOZs, ability to do static-only vs all KOZs
# point to point check against those KOZs

struct Problem
    waypoints::AbstractMatrix{Float64}
    centralKOZDist::Float64
    coneKOZLOS::AbstractMatrix{Float64}
    coneKOZHalfAngle::AbstractVector{Float64}
    sma::Float64
end

function getRandCone(waypoints::AbstractMatrix)
    while true
        RA = rand() * 2 * pi
        DEC = (pi / 12) * randn()
        vec = sph2cart([RA, DEC, 1])
        halfang = rand() * (pi / 12)
        violates = false
        for waypoint in eachcol(waypoints)
            if angVec(vec, waypoint[1:3]) < halfang
                violates = true
                break
            end
        end
        if !violates
            return (vec, halfang)
        end
    end
end

function getBaseProblem()
    waypoints = zeros(6, 5)
    waypoints[:, 1] = [0.0, 1, 0, 0, 0, 0]
    waypoints[:, 2] = [1.0, 0, 0, -4.4329795100374526e-5, -0.0001236164057455293, 0]
    waypoints[:, 3] = [0.0, -1, 0, -2.8278424012935656e-5, -2.2225911392726972e-5, 0]
    waypoints[:, 4] = [-1.0, 0, 0, 4.4329795100374526e-5, 0.0001236164057455293, 0]
    waypoints[:, 5] = [0.0, 1, 0, 0, 0, 0]

    centralDist = 0.1

    coneKOZLOS = zeros(3, 10)
    coneKOZHA = zeros(10)

    temp = sqrt(2) / 2
    temp2 = 10 * pi / 180
    coneKOZLOS[:, 1] = [temp, temp, 0]
    coneKOZHA[1] = temp2

    coneKOZLOS[:, 2] = [temp, -temp, 0]
    coneKOZHA[2] = temp2

    coneKOZLOS[:, 3] = [-temp, -temp, 0]
    coneKOZHA[3] = temp2

    coneKOZLOS[:, 4] = [-temp, temp, 0]
    coneKOZHA[4] = temp2

    for i in 5:10
        coneKOZLOS[:, i], coneKOZHA[i] = getRandCone(waypoints)
    end

    return Problem(waypoints, centralDist, coneKOZLOS, coneKOZHA, 42164.15405)
end

function isFreeSample(
    prob::Problem,
    sample::AbstractVector;
    checkDynamic::Bool=false,
)
    #Checks if the sample is in free space
    if norm(@view(sample[1:3])) < prob.centralKOZDist
        return false
    end

    if checkDynamic
        for i in 1:length(prob.coneKOZHalfAngle)
            if angVec(sample[1:3], prob.coneKOZLOS[:, i]) < probl.coneKOZHalfAngle[i]
                return false
            end
        end
    end

    return true
end

function isFreePath(
    prob::Problem,
    sample::AbstractVector,
    t::AbstractFloat;
    checkDynamic::Bool=false,
)
    x = copy(sample)
    phi = stm(1.0)
    timeElapsed = 0.0

    while timeElapsed < t
        if !isFreeSample(prob, x, checkDynamic=checkDynamic)
            return false
        end
        x = phi * x
        timeElapsed += 1
    end

    phi = stm(t)

    return isFreeSample(prob, phi * sample, checkDynamic=checkDynamic)

end