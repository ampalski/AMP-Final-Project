using JLD2
#prob = getBaseProblem()
#samplingstruct = SampleFree(prob, 400)
#findAllNeighbors!(samplingstruct, prob, .002, 7200.0)

a = load("n400Connections.jld2")
prob = a["prob"]
samplingstruct = a["samplingstruct"]