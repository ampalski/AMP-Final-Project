abstract type Problem end
struct StaticProblem <: Problem
    waypoints::AbstractMatrix{Float64}
    centralKOZDist::Float64
    coneKOZLOS::AbstractMatrix{Float64}
    coneKOZHalfAngle::AbstractVector{Float64}
    sma::Float64
end

mutable struct MutableProblem <: Problem
    waypoints::AbstractMatrix{Float64}
    centralKOZDist::Float64
    coneKOZLOS::AbstractMatrix{Float64}
    coneKOZHalfAngle::AbstractVector{Float64}
    sma::Float64
end

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

abstract type Solution end
mutable struct RRTSoln <: Solution
    solnPath::AbstractMatrix
    totΔv::Float64
    times::AbstractVector
    Δvs::AbstractMatrix
    valid::Bool
end

mutable struct FMTSoln <: Solution
    solnPath::AbstractMatrix
    totΔv::Float64
    times::AbstractVector
    Δvs::AbstractMatrix
    valid::Bool
end