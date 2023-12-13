function SteeringSoln(t::AbstractFloat, x0::AbstractVector, x1::AbstractVector; sma::AbstractFloat=42164.15405)
    if @views(x0[1:3] == x1[1:3]) || t == 0.0
        return vcat(zeros(3), x1[4:6] - x0[4:6])
    end

    stmt = stm(t; sma=sma)
    stm0 = stm(0.0; sma=sma)
    phiv = hcat(stmt[:, 4:6], stm0[:, 4:6])
    right = x1 - stmt * x0

    try
        temp = phiv \ right
    catch
        display(t)
        display(phiv)
        display(right)
    end

    return (phiv \ right)
end

function minTSteeringLineSearch(
    t::AbstractFloat,
    updateSize::AbstractFloat,
    x0::AbstractVector,
    x1::AbstractVector,
    tbounds::Tuple{AbstractFloat,AbstractFloat};
    sma::AbstractFloat=42164.15405,
)
    alpha = 1.0
    dv1 = SteeringSoln(t, x0, x1; sma=sma)
    baseline = norm(dv1[1:3]) + norm(dv1[4:6])
    dv1 = SteeringSoln(t - alpha * updateSize, x0, x1; sma=sma)
    tempUpdate = norm(dv1[1:3]) + norm(dv1[4:6])
    if tempUpdate > baseline
        alphascale = 0.8
        while tempUpdate > baseline
            alpha *= alphascale
            dv1 = SteeringSoln(t - alpha * updateSize, x0, x1; sma=sma)
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
            dv1 = SteeringSoln(t - alpha * updateSize, x0, x1; sma=sma)
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


function minTSteeringOpt(tbounds::Tuple{AbstractFloat,AbstractFloat}, x0::AbstractVector, x1::AbstractVector; sma::AbstractFloat=42164.15405)
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
        dv1 = SteeringSoln(t, x0, x1; sma=sma)
        unperturbed = norm(dv1[1:3]) + norm(dv1[4:6])
        dv1 = SteeringSoln(t + δt, x0, x1; sma=sma)
        perturbed = norm(dv1[1:3]) + norm(dv1[4:6])
        updateSize = (perturbed - unperturbed) / δt
        updateSize = minTSteeringLineSearch(t, updateSize, x0, x1, tbounds; sma=sma)
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

#Find the vector of Δvs to get from root to the specified node
function findCombinedMnvs(problem::Problem, samplingstruct::SamplingMembers, node::Int)
    dvs = Vector{Vector{Float64}}()
    curNode = node
    if isempty(inneighbors(samplingstruct.liveGraph, curNode))
        return dvs
    end

    childNode = curNode
    curNode = inneighbors(samplingstruct.liveGraph, curNode)[1]

    while !isempty(inneighbors(samplingstruct.liveGraph, curNode))
        parentNode = inneighbors(samplingstruct.liveGraph, curNode)[1]
        postMnvState = samplingstruct.samples[curNode] +
                       vcat(zeros(3), samplingstruct.fullEdgeCosts[(curNode, childNode)][1])
        preMnvState = samplingstruct.samples[parentNode] +
                      vcat(zeros(3), samplingstruct.fullEdgeCosts[(parentNode, curNode)][1])
        preMnvState = HCW_Prop(preMnvState, samplingstruct.fullEdgeCosts[(parentNode, curNode)][3], problem.sma)

        pushfirst!(dvs, @views(postMnvState[4:6] - preMnvState[4:6]))

        childNode = curNode
        curNode = parentNode
    end

    pushfirst!(dvs, samplingstruct.fullEdgeCosts[(curNode, childNode)][1])
    return dvs
end

function recalcDV!(soln, prob::Problem)
    total = 0.0
    state = soln.solnPath[:, 1]
    for i in 1:(size(soln.solnPath, 2)-1)
        nextState = soln.solnPath[:, i+1]
        dt = soln.times[i]
        dv = SteeringSoln(dt, state, nextState, sma=prob.sma)[1:3]
        state = state + vcat(zeros(3), dv)
        state = HCW_Prop(state, dt, prob.sma)

        soln.Δvs[:, i] = dv
        total += norm(dv)
    end

    soln.totΔv = total
end

# Given the starting state, the problem waypoints, and the times of maneuvers,
# find the fuel optimal maneuvers, not considering constraints 
function minFuelPath(soln::Solution, prob::Problem)

    tol = 1e-8
    δΔv = 0.0001
    allMnvs = copy(soln.Δvs)
    problemWaypointInd = 2
    curSolnInd = 1
    dvStartInd = 1

    #solve the optimal path segment by segment
    while true
        startWP = soln.solnPath[:, curSolnInd]
        state = Vector{Float64}()
        while @views(soln.solnPath[:, curSolnInd] != prob.waypoints[:, problemWaypointInd])
            state = [state, soln.Δvs[:, curSolnInd]]
            curSolnInd += 1
        end

        #Shooting method
        ktr = 0
        updateSize = 10000000000.0
        while abs(updateSize) > tol && ktr < 1000
            ktr += 1


        end

        problemWaypointInd += 1
        if problemWaypointInd > size(prob.waypoints, 2)
            break
        end
    end

end
#Need to add a flyout utils function, takes a starting state, times of burns and dvs, propagates out to end
#then an error fn that does the flyout and computes an error vector of [total delta v, 
#actual solution probably needs to analytically solve final dv to make sure waypoint is hit exactly