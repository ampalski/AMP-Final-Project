#Need:
# DONE: struct that outlines what configuration the target is in, where the cones are, etc.
# DONE: instantaneous point check against those KOZs, ability to do static-only vs all KOZs
# point to point check against those KOZs

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

function getBaseProblem(; staticProb::Bool=true)
    waypoints = zeros(6, 5)
    waypoints[:, 1] = [0.0, 1, 0, 0, 0, 0]
    waypoints[:, 2] = [1.0, 0, 0, -4.4329795100374526e-5, -0.0001236164057455293, 0]
    waypoints[:, 3] = [0.0, -1, 0, -2.8278424012935656e-5, -2.2225911392726972e-5, 0]
    waypoints[:, 4] = [-1.0, 0, 0, 4.4329795100374526e-5, 0.0001236164057455293, 0]
    waypoints[:, 5] = [0.0, 1, 0, 0, 0, 0]
    # waypoints = zeros(6, 2)
    # waypoints[:, 1] = [0.0, 1, 0, 0, 0, 0]
    # waypoints[:, 2] = [1.0, 0, 0, -4.4329795100374526e-5, -0.0001236164057455293, 0]
    # waypoints[:, 3] = [0.0, -1, 0, -2.8278424012935656e-5, -2.2225911392726972e-5, 0]

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

    if staticProb
        return StaticProblem(waypoints, centralDist, coneKOZLOS, coneKOZHA, 42164.15405)
    else
        return MutableProblem(waypoints, centralDist, coneKOZLOS, coneKOZHA, 42164.15405)
    end
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
            if angVec(sample[1:3], prob.coneKOZLOS[:, i]) < prob.coneKOZHalfAngle[i]
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
    Δt = 0.01 * t
    phi = stm(Δt)
    timeElapsed = 0.0

    while timeElapsed < t
        if !isFreeSample(prob, x, checkDynamic=checkDynamic)
            return false
        end
        x = phi * x
        timeElapsed += Δt
    end

    phi = stm(t)

    return isFreeSample(prob, phi * sample, checkDynamic=checkDynamic)

end

function randRotateProblem!(prob::MutableProblem)
    # Rotate about a random primary axis by a small random amount

    axis = rand(1:3)
    angle = pi / 1800 * randn()
    dcm = zeros(3, 3)
    c = cos(angle)
    s = sin(angle)
    if axis == 1
        dcm[1, 1] = 1.0
        dcm[2, 2] = c
        dcm[3, 3] = c
        dcm[2, 3] = -s
        dcm[3, 2] = s
    elseif axis == 2
        dcm[2, 2] = 1.0
        dcm[1, 1] = c
        dcm[3, 3] = c
        dcm[3, 1] = -s
        dcm[1, 3] = s
    else
        dcm[3, 3] = 1.0
        dcm[2, 2] = c
        dcm[1, 1] = c
        dcm[1, 2] = -s
        dcm[2, 1] = s
    end

    ind = 1
    for vec in eachcol(prob.coneKOZLOS)
        prob.coneKOZLOS[:, ind] = dcm * vec
        ind += 1
    end
end