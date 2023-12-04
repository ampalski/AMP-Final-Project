using GLMakie
crossMatrix(x, y, z) = [0 -z y; z 0 -x; -y x 0]
function plotRRT(
    prob::Problem,
    samplingstruct::SamplingMembers,
    soln::RRTSoln,
)
    set_theme!(theme_black())
    fig = Figure(resolution=(1600, 800))

    ax = Axis3(fig[1, 1])
    ax.xlabel = "Radial"
    ax.ylabel = "In-Track"
    ax.zlabel = "Cross-Track"
    # Plot the KOZs
    waypoints = scatter!(ax, prob.waypoints[1, :], prob.waypoints[2, :], prob.waypoints[3, :], color=:blue, label="Waypoints")
    centralSphere = wireframe!(ax, Sphere(Point3(0.0), prob.centralKOZDist); color=(:white, 0.2))

    n1 = [0, 0, 1.0]
    for coneInd in 1:length(prob.coneKOZHalfAngle)
        dir = prob.coneKOZLOS[:, coneInd]
        halfAng = prob.coneKOZHalfAngle[coneInd]
        dir ./= norm(dir)
        n2 = cross(n1, dir)
        n2 ./= norm(n2)
        rot = acos(dot(dir, n1))
        R = cos(rot) * I(3) + sin(rot) * crossMatrix(n2[1], n2[2], n2[3]) + (1 - cos(rot)) * (n2 * n2')
        radius = tan(halfAng)

        h = 1.0
        u = LinRange(0, h, 25)
        θ = LinRange(0, 2 * pi, 25)
        x = [u * radius * cos(θ) for u in u, θ in θ]
        y = [u * radius * sin(θ) for u in u, θ in θ]
        z = [u for u in u, θ in θ]

        for i in 1:25
            for j in 1:25
                temp = [x[i, j], y[i, j], z[i, j]]
                temp = R * temp
                x[i, j] = temp[1]
                y[i, j] = temp[2]
                z[i, j] = temp[3]
            end
        end

        wireframe!(ax, x, y, z; color=(:white, 0.2))
        #surface!(ax, x, y, z; colormap=:bone_1)
    end

    #Plot the RRT search tree
    helper = zeros(3, 2)
    for i in 1:maximum(keys(samplingstruct.samples))
        helper[:, 1] = samplingstruct.samples[i][1:3]
        for j in outneighbors(samplingstruct.liveGraph, i)
            helper[:, 2] = samplingstruct.samples[j][1:3]
            lines!(ax, helper[1, :], helper[2, :], helper[3, :]; color=:gray)
        end

    end

    # Plot the solution
    if soln.valid
        lines!(ax, soln.solnPath[1, :], soln.solnPath[2, :], soln.solnPath[3, :]; color=:green, linewidth=3.0)
    end

    xlims!(ax, -2.0, 2.0)
    ylims!(ax, -2.0, 2.0)
    zlims!(ax, -2.0, 2.0)
    fig

end


function plotFMT(
    prob::Problem,
    samplingstruct::SamplingMembers,
    soln::FMTSoln,
)
    set_theme!(theme_black())
    fig = Figure(resolution=(1600, 800))

    ax = Axis3(fig[1, 1])
    ax.xlabel = "Radial"
    ax.ylabel = "In-Track"
    ax.zlabel = "Cross-Track"
    # Plot the KOZs
    waypoints = scatter!(ax, prob.waypoints[1, :], prob.waypoints[2, :], prob.waypoints[3, :], color=:blue, label="Waypoints")
    centralSphere = wireframe!(ax, Sphere(Point3(0.0), prob.centralKOZDist); color=(:white, 0.2))

    n1 = [0, 0, 1.0]
    for coneInd in 1:length(prob.coneKOZHalfAngle)
        dir = prob.coneKOZLOS[:, coneInd]
        halfAng = prob.coneKOZHalfAngle[coneInd]
        dir ./= norm(dir)
        n2 = cross(n1, dir)
        n2 ./= norm(n2)
        rot = acos(dot(dir, n1))
        R = cos(rot) * I(3) + sin(rot) * crossMatrix(n2[1], n2[2], n2[3]) + (1 - cos(rot)) * (n2 * n2')
        radius = tan(halfAng)

        h = 1.5
        u = LinRange(0, h, 25)
        θ = LinRange(0, 2 * pi, 25)
        x = [u * radius * cos(θ) for u in u, θ in θ]
        y = [u * radius * sin(θ) for u in u, θ in θ]
        z = [u for u in u, θ in θ]

        for i in 1:25
            for j in 1:25
                temp = [x[i, j], y[i, j], z[i, j]]
                temp = R * temp
                x[i, j] = temp[1]
                y[i, j] = temp[2]
                z[i, j] = temp[3]
            end
        end

        wireframe!(ax, x, y, z; color=(:white, 0.2))
        #surface!(ax, x, y, z; colormap=:bone_1)
    end

    #Plot the FMT search tree
    helper = zeros(3, 10)
    for i in 1:maximum(keys(samplingstruct.samples))
        state = samplingstruct.samples[i]
        helper[:, 1] = state[1:3]
        for j in outneighbors(samplingstruct.liveGraph, i)
            state = samplingstruct.samples[i]
            helper2 = samplingstruct.fullEdgeCosts[(i, j)]
            state += vcat(zeros(3), helper2[1])
            Φ = stm(helper2[3] / 9, sma=prob.sma)
            for k in 2:9
                state = Φ * state
                helper[:, k] = state[1:3]
            end
            helper[:, 10] = samplingstruct.samples[j][1:3]
            lines!(ax, helper[1, :], helper[2, :], helper[3, :]; color=:gray)
        end

    end

    # Plot the solution
    if soln.valid
        scatter!(ax, soln.solnPath[1, :], soln.solnPath[2, :], soln.solnPath[3, :]; color=:green)
        helper = zeros(3, 10)
        state = zeros(6)
        for ind in axes(soln.solnPath, 2)
            if ind == 1
                state = soln.solnPath[:, 1] + vcat(zeros(3), soln.Δvs[:, 1])
            elseif ind == size(soln.solnPath, 2)
                continue
            else
                state += vcat(zeros(3), soln.Δvs[:, ind])
            end
            helper[:, 1] = state[1:3]
            Φ = stm(soln.times[ind] / 9; sma=prob.sma)
            for k in 2:10
                state = Φ * state
                helper[:, k] = state[1:3]
            end

            lines!(ax, helper[1, :], helper[2, :], helper[3, :]; color=:green, linewidth=3.0)
        end
    end

    xlims!(ax, -2.0, 2.0)
    ylims!(ax, -2.0, 2.0)
    zlims!(ax, -2.0, 2.0)
    fig

end