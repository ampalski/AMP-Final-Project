
a = load("RRTNoPost.jld2")
b = load("RRTWiPost.jld2")
c = load("FMTNoPost.jld2")
d = load("FMTWiPost.jld2")

y1 = (1 / 3600) .* c["tof"]
y2 = (1 / 3600) .* d["tof"]
y1 = 1000 .* a["Δvs"]
y2 = 1000 .* c["Δvs"]

y1 = a["times"]
y2 = c["times"]
x1 = fill(1, length(y1))
x2 = fill(2, length(y2))

set_theme!(theme_black())
fig = Figure(resolution=(600, 600))

ax = Axis(fig[1, 1], xticks=([1, 2], ["RRT", "FMT*"]), ylabel="Run Times (s)")

boxplot!(ax, x1, y1)
boxplot!(ax, x2, y2)

fig
############################################################
# a = load("RRTReRun600.jld2")
# b = load("RRTReRun1000.jld2")
# c = load("RRTReRun3600.jld2")
a = load("FMTReRun600.jld2")
b = load("FMTReRun1000.jld2")
c = load("FMTReRun3600.jld2")

y1 = (1 / 3600) .* a["tstep"]
y2 = (1 / 3600) .* b["tstep"]
y3 = (1 / 3600) .* c["tstep"]

x1 = fill(1, length(y1))
x2 = fill(2, length(y2))
x3 = fill(3, length(y3))

set_theme!(theme_black())
fig = Figure(resolution=(600, 600))

ax = Axis(fig[1, 1], xticks=([1, 2, 3], ["600", "1000", "3600"]), xlabel="Execution Rate (s)",
    ylabel="Time Until KOZ Violated or Goal Met (hrs)")

boxplot!(ax, x1, y1)
boxplot!(ax, x2, y2)
boxplot!(ax, x3, y3)

ylims!(ax, [0, 24])

fig