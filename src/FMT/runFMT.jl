using JLD2

function runFMT(numRuns::Int)
    times = zeros(numRuns)
    ΔvVec = zeros(numRuns)
    successes = Vector{Bool}()
    tof = zeros(numRuns)

    for i in 1:numRuns
        display(i)

        a = load("n2000Connections.jld2")
        prob = getBaseProblem()
        samplingstruct = a["samplingstruct"]
        soln = getEmptyFMTSoln()
        tStart = time()
        FMTPlan!(samplingstruct, soln, prob)
        times[i] = time() - tStart
        if soln.valid
            push!(successes, true)
            ΔvVec[i] = soln.totΔv
            tof[i] = sum(soln.times)
        else
            push!(successes, false)
        end
    end

    return (times, ΔvVec, successes, tof)
end

times, Δvs, successes, tof = runFMT(100)