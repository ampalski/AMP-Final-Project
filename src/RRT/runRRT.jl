#= 
function that runs RRT n times
record the number of successes
runTime
delta-v 
# of waypoints?
time of arrival?
=#


function runRRT(numRuns::Int)
    times = zeros(numRuns)
    ΔvVec = zeros(numRuns)
    successes = Vector{Bool}()
    tof = zeros(numRuns)

    for i in 1:numRuns
        display(i)
        prob = getBaseProblem()
        samplingstruct = getEmptySampleStruct(prob, 2000)
        soln = getEmptyRRTSoln()
        tStart = time()
        RRTPlan!(samplingstruct, soln, prob)
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

times, Δvs, successes, tof = runRRT(100)