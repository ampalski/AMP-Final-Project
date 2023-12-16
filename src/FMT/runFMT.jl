using JLD2

function runFMT(numRuns::Int; postProcess::Bool=false)
    times = zeros(numRuns)
    ΔvVec = zeros(numRuns)
    successes = Vector{Bool}()
    tof = zeros(numRuns)
    t_waypoints = [21600.0, 21600.0, 21600.0, 21600.0]

    for i in 1:numRuns
        display(i)

        a = load("n2000Connections.jld2")
        prob = getBaseProblem()
        samplingstruct = a["samplingstruct"]
        soln = getEmptyFMTSoln()
        tStart = time()
        FMTPlan!(samplingstruct, soln, prob)
        display("Planning Done")

        if soln.valid
            if postProcess
                getCorrectTimes!(soln, prob, t_waypoints)
                recalcDV!(soln, prob)
                display("Time Correction Done")
            end
            times[i] = time() - tStart
            push!(successes, true)
            ΔvVec[i] = soln.totΔv
            tof[i] = sum(soln.times)
        else
            times[i] = time() - tStart
            push!(successes, false)
        end
    end

    return (times, ΔvVec, successes, tof)
end

# times, Δvs, successes, tof = runFMT(100)
times, Δvs, successes, tof = runFMT(100; postProcess=true)


function runRealTimeFMT(runRate::AbstractFloat)
    logfile = FormatLogger(open("log.txt", "w")) do io, args
        # Write the module, level and message only
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end

    prob = getBaseProblem(staticProb=false)
    a = load("n2000Connections.jld2")
    samplingstruct = a["samplingstruct"]
    rerun = true

    state = prob.waypoints[:, 1]
    dt = 1.0
    curdt = dt
    Φ = stm(dt; sma=prob.sma)
    curΦ = copy(Φ)
    t = 0.0
    t_waypoints = [21600.0, 21600.0, 21600.0, 21600.0]
    tEnd = 1.5 * sum(t_waypoints)
    success = true
    path = zeros(4, 2 * Int(sum(t_waypoints)))
    pathInd = 1
    tToRerun = runRate
    tleft = runRate
    soln = getEmptyFMTSoln()
    onStep = 1
    totaldv = 0.0

    while success
        #Check to see if new pathfinding is required, run FMT if soln
        if abs(t % 10) < 0.999
            with_logger(logfile) do
                #display(t)
                @debug t
            end
        end
        if rerun && t_waypoints[1] > runRate
            with_logger(logfile) do
                #display(t)
                @debug "Running FMT*"
                @debug @sprintf("Time to next problem waypoint: %f", t_waypoints[1])
            end
            prob.waypoints[:, 1] = state
            addOnlineNode!(samplingstruct, state, 0.001, tleft)

            soln = getEmptyFMTSoln()
            ktr2 = 0
            while !soln.valid && ktr2 < 20
                ktr2 += 1
                n = nv(samplingstruct.liveGraph)
                samplingstruct.liveGraph = SimpleDiGraph(n)
                empty!(samplingstruct.unvisitedVertex)
                for i in 1:n
                    i == samplingstruct.nodeInit && continue
                    push!(samplingstruct.unvisitedVertex, i)
                end
                empty!(samplingstruct.openVertex)
                push!(samplingstruct.openVertex, samplingstruct.nodeInit)
                empty!(samplingstruct.closedVertex)
                FMTPlan!(samplingstruct, soln, prob)
            end
            if ktr2 >= 20
                success = false
                with_logger(logfile) do
                    @debug "FMT* failed, failing out."
                end
                continue
            end
            ktr2 = 0
            while abs(sum(soln.times) - sum(t_waypoints)) > 1.0 && ktr2 < 20
                ktr2 += 1
                getCorrectTimes!(soln, prob, t_waypoints)
            end
            recalcDV!(soln, prob)
            onStep = 2
            tleft = soln.times[1]
            tToRerun = runRate
            state += vcat(zeros(3), soln.Δvs[:, 1])
            totaldv += norm(soln.Δvs[:, 1])

            with_logger(logfile) do
                #display(t)
                @debug @sprintf("Time to next solution waypoint: %f", tleft)
                @debug @sprintf("t_waypoints time left: %f", sum(t_waypoints))
                @debug @sprintf("solution time left: %f", sum(soln.times))
            end

            rerun = false
        end

        #Step forward in time
        state = curΦ * state
        t += curdt
        tToRerun -= curdt
        tleft -= curdt
        t_waypoints[1] -= curdt
        if t_waypoints[1] < 0
            with_logger(logfile) do
                #display(t)
                @debug @sprintf("t=%f: %f, %f, %f", t, state[1], state[2], state[3])
            end
        end

        path[1, pathInd] = t
        path[2:4, pathInd] = state[1:3]
        pathInd += 1

        #Check for waypoints
        if norm(state[1:3] - soln.solnPath[1:3, onStep]) < 0.00001 || abs(tleft) < 0.001
            with_logger(logfile) do
                #display(t)
                @debug @sprintf("Hit solution waypoint at t=%f", t)
            end

            if onStep <= size(soln.Δvs, 2)
                state += vcat(zeros(3), soln.Δvs[:, onStep])
                totaldv += norm(soln.Δvs[:, onStep])
                tleft = soln.times[onStep]
                onStep += 1
            else
                tleft = runRate
                success = false
                with_logger(logfile) do
                    #display(t)
                    @debug "Hit solution but no more delta-vs"
                    @debug state
                end
                display(state)
            end

            with_logger(logfile) do
                #display(t)
                @debug @sprintf("Time to next solution waypoint: %f", tleft)
            end
        end

        if norm(state[1:3] - prob.waypoints[1:3, 2]) < 0.0001
            if size(prob.waypoints, 2) == 2
                break
            end
            with_logger(logfile) do
                #display(t)
                @debug @sprintf("Hit problem waypoint at t=%f", t)
            end
            prob.waypoints = prob.waypoints[:, 2:end]
            t_waypoints = t_waypoints[2:end]
            samplingstruct.nodeInit += 1
        end

        #Check for timestep handling

        if tleft < dt
            with_logger(logfile) do
                @debug @sprintf("Resetting curDT: %f", tleft)
            end

            curdt = tleft
            curΦ = stm(curdt; sma=prob.sma)
        elseif curdt < dt
            with_logger(logfile) do
                @debug "Resetting curDT back to normal"
            end
            curdt = dt
            curΦ = copy(Φ)
        end

        if tToRerun < 0
            rerun = true
            tToRerun = runRate
        end

        # Rotate problem until waypoints are not covered
        isValid = false
        while !isValid
            randRotateProblem!(prob)
            isValid = true
            for waypoint in eachcol(prob.waypoints)
                for i in 1:length(prob.coneKOZHalfAngle)
                    if angVec(waypoint[1:3], prob.coneKOZLOS[:, i]) < prob.coneKOZHalfAngle[i]
                        isValid = false
                        break
                    end
                end
                if !isValid
                    break
                end
            end
        end

        if !isFreeSample(prob, state; checkDynamic=true) || t > tEnd
            success = false
            with_logger(logfile) do
                @debug "Violated KOZ, failing out."
            end
        end

    end
    return (success, path[:, 1:pathInd], totaldv)
end

success = Vector{Bool}()
tstep = zeros(10)
dvs = zeros(10)
for i in 1:10
    temp, path, dvs[i] = runRealTimeFMT(600.0)
    if norm(path[2:4, end-1] - [0.0, 1, 0]) < 0.01
        temp = true
    end
    push!(success, temp)
    tstep[i] = path[1, end-1]
end
dvs .*= 1000.0