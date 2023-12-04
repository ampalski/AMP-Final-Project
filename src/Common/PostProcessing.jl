# Things to add: 
# Selecting an edge at random and attempting to change the time to get to tf
# Finding optimal, unconstrained n-burn sequence
# Morphing trajectory towards unconstrained to min delta-v

function getCorrectTimes!(soln, prob::Problem, t_waypoints::AbstractVector)
    #Find the points in soln.solnPath that correspond to the problem waypoints
    #For each problem waypoint, adjust the intermediate points to align 
    #the time to the waypoint to the corresponding value in t_waypoints

    #Take the current difference between t_waypoints and current time-to-waypoint
    #Take a random intermediate path
    #randomly choose between:
    #Propagate for some fraction of the difference more/less (moving the waypoint)
    #Change the time to waypoint by some fraction of the difference
    #Until difference goes away
    #Then go to next waypoint
end