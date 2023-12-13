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
    prob = getBaseProblem(staticProb=false)

end