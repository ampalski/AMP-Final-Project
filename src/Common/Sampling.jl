# At least need a SampleFree(n) -- have sample free only check against static KOZs
# Probably wouldn't hurt to have something that samples in a gaussian about a specific point with a specific covariance

mutable struct SamplingMembers
    feasGraph::Graphs.SimpleGraphs.AbstractSimpleGraph{Int64}
    liveGraph::Graphs.SimpleGraphs.AbstractSimpleGraph{Int64}
    edgeCosts::SparseArrays.AbstractSparseMatrixCSC{Float64,Int64}
    fullEdgeCosts::Dict{Tuple{Int,Int},Tuple{Vector{Float64},Vector{Float64},Float64}}
    samples::Dict{Int,Vector{Float64}}
    openVertex::Set{Int} #Try these as sets, may need to change if re-ordering later
    unvisitedVertex::Set{Int} #Maybe linkedlists.jl, if adding/subtracting is more important than checking if a member is present
    closedVertex::Set{Int}
    nodeInit::Int
end

function SampleFree(prob::Problem, n::Int)
    HaltonSamples = HaltonPoint(6, length=n)
    posBounds = (-2.0, 2.0)
    velBounds = (-0.005, 0.005)
    scale = vcat(repeat([posBounds[2] - posBounds[1]], 3), repeat([velBounds[2] - velBounds[1]], 3))
    shift = vcat(repeat([posBounds[1]], 3), repeat([velBounds[1]], 3))

    samples = Dict{Int,Vector{Float64}}()
    sizehint!(samples, n + 5)
    openVertex = Set{Int}()
    unvisitedVertex = Set{Int}()
    closedVertex = Set{Int}()
    i = 1
    for s in HaltonSamples
        samples[i] = scale .* s .+ shift
        push!(unvisitedVertex, i)
        i += 1
    end

    nodeInit = i
    for s in eachcol(prob.waypoints)
        samples[i] = s
        if i == nodeInit
            push!(openVertex, i)
        else
            push!(unvisitedVertex, i)
        end
        i += 1
    end

    graph = SimpleDiGraph(i - 1)
    liveGraph = SimpleDiGraph(i - 1)

    edgeCosts = spzeros(i - 1, i - 1)

    fullEdgeCosts = Dict{Tuple{Int,Int},Tuple{Vector{Float64},Vector{Float64},Float64}}()

    return SamplingMembers(graph, liveGraph, edgeCosts, fullEdgeCosts, samples, openVertex, unvisitedVertex, closedVertex, nodeInit)
end

function findAllNeighbors!(
    samplingHelper::SamplingMembers,
    prob::Problem,
    dvMax::AbstractFloat,
    tMax::AbstractFloat,
)
    n = length(samplingHelper.samples)
    tbounds = (0.0, tMax)
    for i in 1:n
        i % 100 == 0 && display(i)
        for j in 1:n
            if i == j
                continue
            end
            t = minTSteeringOpt(tbounds, samplingHelper.samples[i], samplingHelper.samples[j])
            dvs = SteeringSoln(t, samplingHelper.samples[i], samplingHelper.samples[j])
            dv = @views(norm(dvs[1:3]) + norm(dvs[4:6]))

            if dv > dvMax
                continue
            end

            x0 = samplingHelper.samples[i] + vcat(zeros(3), dvs[1:3])
            if !isFreePath(prob, x0, t, checkDynamic=false)
                continue
            end

            add_edge!(samplingHelper.feasGraph, i, j)
            samplingHelper.edgeCosts[i, j] = dv
            samplingHelper.fullEdgeCosts[(i, j)] = (dvs[1:3], dvs[4:6], t)

        end
    end

end

