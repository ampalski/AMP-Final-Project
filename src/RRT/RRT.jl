#=
1) Init tree with root = q_init
while q_goal not found:
2) generate random sample
    a) if rand() < p, q_rand = q_goal
    b) else, collisionfreesample
3) find closest node to q_rand, designate q_near
4) generate path from q_near to q_rand
5) take step along that path according to step size r
    a) if collision free, add as q_new
    b) continue stepping until q_rand is reached or collision is found
6) add all collision free paths to tree
7) If (q_new-q_goal).norm() < epsilon, consider goal reached and return
=#

#cpp version has a "findclosestnode" which should already be DONE
#and a "plan" function that did most of the work
#may need some kind of solution struct.

mutable struct RRTSoln
    solnPath::AbstractMatrix
    totΔv::Float64
    times::AbstractVector
    Δvs::AbstractMatrix
    valid::Bool
end

function RRTPlan!(sampleStruct::SamplingMembers, soln::RRTSoln, problem::Problem)
    solnFound = false
    stepDone = false
    currentGoal::Int = 2
    totalGoal::Int = size(problem.waypoints, 2)
    baseNode::Int = sampleStruct.nodeInit
    currentNode::Int = baseNode + 1
    goalNode::Int = 0
    tbounds = (0.0, 21600.0)
    ΔvMax = 0.01
    maxNodes = nv(sampleStruct.liveGraph)

    while !solnFound
        # Generate Random sample
        q_rand = rand() < 0.1 ? problem.waypoints[:, currentGoal] : getRandomSample(problem)

        # Find closest node
        q_near_ind = findClosestNode(q_rand, sampleStruct.samples, minNode=baseNode)
        q_near = sampleStruct.samples[q_near_ind]

        # Create path from existing node to random sample
        t = minTSteeringOpt(tbounds, q_near, q_rand)
        dvs = SteeringSoln(t, q_near, q_rand)
        dv = @views(norm(dvs[1:3]))

        if dv > ΔvMax
            continue
        end

        stepDone = false
        curState = q_near + vcat(zeros(3), dvs[1:3])
        incT = t / 10
        Φ = stm(incT; sma=problem.sma)
        curT = t / 10
        # Step towards q_rand, adding to the tree as it goes
        while !stepDone
            abs(curT - t) < t / 20 && (stepDone = true)
            currentNode == maxNodes && (stepDone = true)

            if !isFreePath(problem, curState, incT; checkDynamic=true)
                stepDone = true
                continue
            end

            newState = Φ * curState
            sampleStruct.samples[currentNode] = newState
            add_edge!(sampleStruct.liveGraph, q_near_ind, currentNode)
            if curT == incT
                sampleStruct.edgeCosts[q_near_ind, currentNode] = dv
                sampleStruct.fullEdgeCosts[(q_near_ind, currentNode)] = (dvs[1:3], zeros(3), incT)
            else
                sampleStruct.edgeCosts[q_near_ind, currentNode] = 0
                sampleStruct.fullEdgeCosts[(q_near_ind, currentNode)] = (zeros(3), zeros(3), incT)
            end
            q_near_ind = currentNode
            curState = copy(newState)
            currentNode += 1
            curT += incT

            # Check for goal
            if @views(norm(newState[1:3] - problem.waypoints[1:3, currentGoal])) < 0.01
                currentGoal += 1
                baseNode = currentNode - 1
                if currentGoal > totalGoal
                    solnFound = true
                    goalNode = baseNode
                end
                break
            end
        end

        if currentNode > maxNodes
            break
        end

    end
    if !solnFound
        soln.valid = false
        return
    end

    # Convert solution
    nodePath = [goalNode]
    currentNode = goalNode
    while currentNode != 1
        currentNode = inneighbors(sampleStruct.liveGraph, currentNode)[1]
        pushfirst!(nodePath, currentNode)
    end
    n = length(nodePath)
    solnPath = zeros(6, n)
    times = zeros(n - 1)
    Δvs = zeros(3, n - 1)
    totΔv = 0.0
    for i in 1:n
        solnPath[:, i] = sampleStruct.samples[nodePath[i]]
        if i == n
            break
        end

        allEdgeCosts = sampleStruct.fullEdgeCosts[(nodePath[i], nodePath[i+1])]
        times[i] = allEdgeCosts[3]
        Δvs[:, i] = allEdgeCosts[1]
        totΔv += norm(allEdgeCosts[1])
    end

    soln.valid = true
    soln.solnPath = solnPath
    soln.totΔv = totΔv
    soln.times = times
    soln.Δvs = Δvs
    return
end

function getEmptyRRTSoln()
    return RRTSoln(zeros(2, 2), 0.0, zeros(2), zeros(2, 2), false)
end