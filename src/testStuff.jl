timer = time()
prob = getBaseProblem()
samplingstruct = SampleFree(prob, 2000)
findAllNeighbors!(samplingstruct, prob, 0.001, 10800.0)
timer = time() - timer

using JLD2
a = load("n2000Connections.jld2")
prob = getBaseProblem()
samplingstruct = a["samplingstruct"]
soln = getEmptyFMTSoln()
FMTPlan!(samplingstruct, soln, prob)
display(soln.valid)
# plotFMT(prob, samplingstruct, soln)
soln2 = deepcopy(soln)
getCorrectTimes!(soln, prob, [21600.0, 21600.0, 21600.0, 21600.0])
recalcDV!(soln, prob)
plotPostProcessed(prob, soln2, soln)


samplingstruct = a["samplingstruct"]
soln = getEmptyFMTSoln()
@run FMTPlan!(samplingstruct, soln, prob)

prob = getBaseProblem()
samplingstruct = getEmptySampleStruct(prob, 2000)
soln = getEmptyRRTSoln()
while !soln.valid
    samplingstruct = getEmptySampleStruct(prob, 2000)
    RRTPlan!(samplingstruct, soln, prob)
end
display(soln.valid)
length(samplingstruct.samples)
soln2 = deepcopy(soln)
getCorrectTimes!(soln, prob, [21600.0, 21600.0, 21600.0, 21600.0])
recalcDV!(soln, prob)
plotPostProcessed(prob, soln2, soln)
plotRRT(prob, samplingstruct, soln)


logfile = FormatLogger(open("log.txt", "w")) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
temp = 123
with_logger(logfile) do
    #display(t)
    @debug @sprintf("Testing %d", temp)
end


times, Δvs, successes, tof = runRRT(100)
jldsave("RRTNoPost.jld2"; times, Δvs, successes, tof)
times, Δvs, successes, tof = runRRT(100; postProcess=true)
jldsave("RRTWiPost.jld2"; times, Δvs, successes, tof)

times, Δvs, successes, tof = runFMT(100)
jldsave("FMTNoPost.jld2"; times, Δvs, successes, tof)
times, Δvs, successes, tof = runFMT(100; postProcess=true)
jldsave("FMTWiPost.jld2"; times, Δvs, successes, tof)