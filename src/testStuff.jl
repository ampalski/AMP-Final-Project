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
plotFMT(prob, samplingstruct, soln)

samplingstruct = a["samplingstruct"]
soln = getEmptyFMTSoln()
@run FMTPlan!(samplingstruct, soln, prob)

prob = getBaseProblem()
samplingstruct = getEmptySampleStruct(prob, 2000)
soln = getEmptyRRTSoln()
RRTPlan!(samplingstruct, soln, prob)

display(soln.valid)
length(samplingstruct.samples)

plotRRT(prob, samplingstruct, soln)