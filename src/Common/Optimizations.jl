function SteeringSoln(t::AbstractFloat, x0::AbstractVector, x1::AbstractVector)
    if @views(x0[1:3] == x1[1:3])
        return vcat(zeros(3), x1[4:6] - x0[4:6])
    end

    stmt = stm(t)
    stm0 = stm(0.0)
    phiv = hcat(stmt[:, 4:6], stm0[:, 4:6])
    right = x1 - stmt * x0
    return (phiv \ right)
end

function minTSteeringLineSearch(
    t::AbstractFloat,
    updateSize::AbstractFloat,
    x0::AbstractVector,
    x1::AbstractVector,
    tbounds::Tuple{AbstractFloat,AbstractFloat}
)
    alpha = 1.0
    dv1 = SteeringSoln(t, x0, x1)
    baseline = norm(dv1[1:3]) + norm(dv1[4:6])
    dv1 = SteeringSoln(t - alpha * updateSize, x0, x1)
    tempUpdate = norm(dv1[1:3]) + norm(dv1[4:6])
    if tempUpdate > baseline
        alphascale = 0.8
        while tempUpdate > baseline
            alpha *= alphascale
            dv1 = SteeringSoln(t - alpha * updateSize, x0, x1)
            tempUpdate = norm(dv1[1:3]) + norm(dv1[4:6])

            if alpha < 0.0001
                break
            end
        end
    else
        alphascale = 1.5
        oldUpdate = tempUpdate
        while true
            alpha *= 1.5
            dv1 = SteeringSoln(t - alpha * updateSize, x0, x1)
            tempUpdate = norm(dv1[1:3]) + norm(dv1[4:6])
            if tempUpdate > oldUpdate
                alpha /= 1.5
                break
            end
            oldUpdate = tempUpdate
            newT = t - alpha * updateSize
            if alpha > 1e10 || newT < tbounds[1] || newT > tbounds[2]
                alpha /= 1.5
                break
            end
        end
    end
    return alpha * updateSize
end


function minTSteeringOpt(tbounds::Tuple{AbstractFloat,AbstractFloat}, x0::AbstractVector, x1::AbstractVector)
    if @views(x0[1:3] == x1[1:3])
        return 0.0
    end

    updateSize = 1000000.0
    tol = 1e-8
    ktr = 0
    t = (tbounds[1] + tbounds[2]) / 2
    δt = 1.0
    while abs(updateSize) > tol
        ktr += 1
        dv1 = SteeringSoln(t, x0, x1)
        unperturbed = norm(dv1[1:3]) + norm(dv1[4:6])
        dv1 = SteeringSoln(t + δt, x0, x1)
        perturbed = norm(dv1[1:3]) + norm(dv1[4:6])
        updateSize = (perturbed - unperturbed) / δt
        updateSize = minTSteeringLineSearch(t, updateSize, x0, x1, tbounds)
        t = t - updateSize
        if t < tbounds[1] || t > tbounds[2]
            t = clamp(t, tbounds[1], tbounds[2])
            break
        end
        if ktr > 1000
            break
        end
    end
    return t
end