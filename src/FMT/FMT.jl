#=
1: Add xinit to the root of the tree T , as a member of the frontier set Vopen
2: Generate samples S ← SampleFree(X , n, t0) and add them to the unexplored set Vunvisited
3: Set the minimum cost-to-come node in the frontier set as z ← xinit
4: while true
5: for each neighbor x of z in Vunvisited
6:      Find the neighbor xmin in Vopen of cheapest cost-to-go to x
7:      Compute the trajectory between them as [x(t), u(t), t] ← Steer(xmin, x) (Section 7.1)
8:      if CollisionFree(x(t), u(t), t)
9:          Add the trajectory from xmin to x to tree T
10: Remove all x from the unexplored set Vunvisited
11: Add any new connections x to the frontier Vopen
12: Remove z from the frontier Vopen and add it to Vclosed
13: if Vopen is empty
14:     return Failure
15: Reassign z as the node in Vopen with smallest cost-to-come from the root (xinit)
16: if z is in the goal region Xgoal
17:     return Success, and the unique trajectory from the root (xinit) to z
=#

function FMTPlan!(samplingstruct::SamplingMembers, soln::FMTSoln, problem::Problem)
    # logfile = FormatLogger(open("log2.txt", "w")) do io, args
    #     # Write the module, level and message only
    #     println(io, args._module, " | ", "[", args.level, "] ", args.message)
    # end

    solnFound = false
    currentGoal = 2
    totalGoal = size(problem.waypoints, 2)
    minCostNode = samplingstruct.nodeInit
    curNodeInit = samplingstruct.nodeInit
    while !solnFound
        #display(length(outneighbors(samplingstruct.feasGraph, minCostNode)))
        # with_logger(logfile) do
        #     @debug @sprintf("Searching outneighbors of %d", minCostNode)
        # end
        for x in outneighbors(samplingstruct.feasGraph, minCostNode)
            if !(x in samplingstruct.unvisitedVertex)
                continue
            end
            minCost = 9999999999.0
            minCostInd = 0
            for xtilde in samplingstruct.openVertex
                if !(xtilde in inneighbors(samplingstruct.feasGraph, x))
                    continue
                end
                cost = samplingstruct.edgeCosts[xtilde, x]
                tempState = samplingstruct.samples[xtilde] + vcat(zeros(3), samplingstruct.fullEdgeCosts[(xtilde, x)][1])
                if cost < minCost && isFreePath(problem, tempState, samplingstruct.fullEdgeCosts[(xtilde, x)][3], checkDynamic=true)
                    minCost = cost
                    minCostInd = xtilde
                end
                # tempState = samplingstruct.samples[xtilde] + vcat(zeros(3), samplingstruct.fullEdgeCosts[(xtilde, x)][1])
                # if isFreePath(problem, tempState, samplingstruct.fullEdgeCosts[(xtilde, x)][3], checkDynamic=true)
                #     add_edge!(samplingstruct.liveGraph, xtilde, x)
                #     push!(samplingstruct.openVertex, x)
                # end
            end

            if minCostInd != 0
                add_edge!(samplingstruct.liveGraph, minCostInd, x)
                push!(samplingstruct.openVertex, x)
            end
            delete!(samplingstruct.unvisitedVertex, x)
        end
        # with_logger(logfile) do
        #     @debug @sprintf("Done")
        # end
        delete!(samplingstruct.openVertex, minCostNode)
        push!(samplingstruct.closedVertex, minCostNode)
        if isempty(samplingstruct.openVertex)
            soln.valid = false
            return
        end
        mindv = 9999999
        mindvInd = 0
        for z in samplingstruct.openVertex
            # with_logger(logfile) do
            #     @debug @sprintf("Finding total dv for %d", z)
            # end
            dvs = findCombinedMnvs(problem, samplingstruct, z)
            # with_logger(logfile) do
            #     @debug @sprintf("Done")
            # end
            dv = sum(norm.(dvs))
            if dv < mindv
                mindv = dv
                mindvInd = z
            end
        end
        minCostNode = mindvInd
        #display("New Min Cost Node: ")
        #display(minCostNode)
        # Check for goal
        if @views(norm(samplingstruct.samples[minCostNode][1:3] - problem.waypoints[1:3, currentGoal])) < 0.01
            # with_logger(logfile) do
            #     @debug @sprintf("Building output for %d", minCostNode)
            # end
            currentGoal += 1
            empty!(samplingstruct.openVertex)
            push!(samplingstruct.openVertex, minCostNode)
            empty!(samplingstruct.closedVertex)

            #Construct output
            nodePath = [minCostNode]
            currentNode = minCostNode
            #display("Constructing Output")
            while currentNode != curNodeInit
                currentNode = inneighbors(samplingstruct.liveGraph, currentNode)[1]
                pushfirst!(nodePath, currentNode)
            end
            #display("Complete")
            totΔv = 0.0
            Δvs = findCombinedMnvs(problem, samplingstruct, minCostNode)
            Δvs2 = zeros(3, length(Δvs))

            for i in axes(Δvs, 1)
                Δvs2[:, i] = Δvs[i]
            end

            n = length(nodePath)
            solnPath = zeros(6, n)
            times = zeros(n - 1)

            for i in 1:n
                solnPath[:, i] = samplingstruct.samples[nodePath[i]]
                if i == n
                    break
                end

                times[i] = samplingstruct.fullEdgeCosts[(nodePath[i], nodePath[i+1])][3]
                totΔv += norm(Δvs[i])
            end
            # with_logger(logfile) do
            #     @debug @sprintf("Saving")
            # end
            # Save output
            if curNodeInit == samplingstruct.nodeInit #first time through
                soln.solnPath = solnPath
                soln.totΔv = totΔv
                soln.times = times
                soln.Δvs = Δvs2
            else
                postMnvState = solnPath[:, 1] + vcat(zeros(3), Δvs2[:, 1])
                temp = SteeringSoln(soln.times[end], soln.solnPath[:, end-1], soln.solnPath[:, end])
                preMnvState = soln.solnPath[:, end] - vcat(zeros(3), temp[4:6])

                totΔv -= norm(Δvs2[:, 1])
                Δvs2[:, 1] = @views(postMnvState[4:6] - preMnvState[4:6])
                totΔv += norm(Δvs2[:, 1])

                soln.solnPath = hcat(soln.solnPath, @view(solnPath[:, 2:end]))
                soln.times = vcat(soln.times, times)
                soln.totΔv += totΔv
                soln.Δvs = hcat(soln.Δvs, Δvs2)
            end
            # with_logger(logfile) do
            #     @debug @sprintf("Resetting graphs")
            # end
            n = nv(samplingstruct.liveGraph)
            samplingstruct.liveGraph = SimpleDiGraph(n)
            for i in 1:n
                if i == minCostNode
                    continue
                end
                push!(samplingstruct.unvisitedVertex, i)
            end
            curNodeInit = minCostNode
            if currentGoal > totalGoal
                solnFound = true
            end
        end
    end
    if solnFound
        soln.valid = true
    end
    return
end


#=
DONE: for each neighbor x of z in Vunvisited
DONE:      Find the neighbor xmin in Vopen of cheapest cost-to-go to x
DONE ish:      Compute the trajectory between them as [x(t), u(t), t] ← Steer(xmin, x) (Section 7.1)
DONE:      if CollisionFree(x(t), u(t), t)
DONE:          Add the trajectory from xmin to x to tree T
DONE: Remove all x from the unexplored set Vunvisited
DONE: Add any new connections x to the frontier Vopen
DONE: Remove z from the frontier Vopen and add it to Vclosed
DONE: if Vopen is empty
DONE:     return Failure
DONE: Reassign z as the node in Vopen with smallest cost-to-come from the root (xinit)
DONE: if z is in the goal region Xgoal
17:     return Success, and the unique trajectory from the root (xinit) to z
=#

function getEmptyFMTSoln()
    return FMTSoln(zeros(2, 2), 0.0, zeros(2), zeros(2, 2), false)
end

