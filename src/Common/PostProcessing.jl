# Things to add: 
# Selecting an edge at random and attempting to change the time to get to tf
# Finding optimal, unconstrained n-burn sequence
# Morphing trajectory towards unconstrained to min delta-v

function getCorrectTimes!(soln, prob::Problem, t_waypoints::AbstractVector)
    #Find the points in soln.solnPath that correspond to the problem waypoints
    #For each problem waypoint, adjust the intermediate points to align 
    #the time to the waypoint to the corresponding value in t_waypoints
    problemWaypointInd = 2
    curSolnInd = 1
    while problemWaypointInd <= size(prob.waypoints, 2)
        inds = [curSolnInd]
        waypoints = soln.solnPath[:, curSolnInd]
        wpTimes = Vector{Float64}()
        while @views(norm(soln.solnPath[1:3, curSolnInd] - prob.waypoints[1:3, problemWaypointInd]) > 0.00001)
            if curSolnInd > (length(soln.times))
                break
            end

            push!(wpTimes, soln.times[curSolnInd])
            curSolnInd += 1
            waypoints = [waypoints soln.solnPath[:, curSolnInd]]
            push!(inds, curSolnInd)
            if curSolnInd > (length(soln.times))
                break
            end
        end

        if isempty(wpTimes)
            return
        end

        Δt = t_waypoints[problemWaypointInd-1] - sum(wpTimes)
        temp = abs(t_waypoints[problemWaypointInd-1]) / length(wpTimes)
        minT = temp > 600 ? 600 : temp

        ktr = 0
        while abs(Δt) > 1.0 && ktr <= 100
            ktr += 1
            changeAmt = rand() < 0.25 ? Δt : Δt * rand()

            if rand() < 0.5
                changeInd = rand(1:length(wpTimes))
            else
                changeInd = changeAmt > 0 ? argmin(wpTimes) : argmax(wpTimes)
            end

            if (wpTimes[changeInd] + changeAmt) < minT
                changeAmt = minT - wpTimes[changeInd]
            end

            if changeAmt == 0.0
                continue
            end

            startWP = waypoints[:, changeInd]
            newTime = changeAmt + wpTimes[changeInd]
            finWP = waypoints[:, changeInd+1]
            if rand() < 0.5 || changeInd == length(wpTimes)
                #Attempt to change the time to existing waypoint    
                dv = SteeringSoln(newTime, startWP, finWP)
                mode = 1
            else
                # Propagate towards original waypoint, with new tof
                dv = SteeringSoln(wpTimes[changeInd], startWP, finWP; sma=prob.sma)
                mode = 2
            end
            startWP = startWP + vcat(zeros(3), dv[1:3])

            if !isFreePath(prob, startWP, newTime; checkDynamic=true)
                continue
            end
            updateInd = inds[changeInd]

            if mode == 2
                finWP = HCW_Prop(startWP, newTime, prob.sma)
                nextWP = waypoints[:, changeInd+2]
                nextTime = wpTimes[changeInd+1]
                dv = SteeringSoln(nextTime, finWP, nextWP)
                finWP = finWP + vcat(zeros(3), dv[1:3])
                if !isFreePath(prob, finWP, nextTime; checkDynamic=true)
                    continue
                end

                waypoints[:, changeInd+1] = finWP
                soln.solnPath[:, updateInd+1] = finWP
            end
            # Update the soln
            wpTimes[changeInd] = newTime
            soln.times[updateInd] = newTime
            Δt -= changeAmt

        end
        problemWaypointInd += 1

    end

end
