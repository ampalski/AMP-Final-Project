General flow:
* Sample n points
* do pre-computation on unbounded connections with limiting delta-v and tof
* Run through the base technique, save off results
* Send results through smoothing process, save off results
* Compare results... accuracy, delta-v, runtimes, etc.

Data structures:
* Graph holds which nodes are connected
* 3x sparse array holds delta-v cost for edge costs
* map? (int, int) -> (vector, vector, float)? to store dv1, dv2, T for connections
* Map holds int -> vector{int} for which other nodes are in its neighborhood after pre-processing
* Map holds sample locations (int, node number -> vector RIC/RICdot)
* Simple array (vector) holds the Vector_open, Vector_unvisited, other sets. Just need to storey node number as int
* Probably have a struct hold all the common features, separate structs for any technique specific data
* Primary waypoint paths can just come back as matrices. But want the ability to get matrices for every path held in the graph for plotting

Things I need to figure out:
* Linear programming. Make sure I can do the min-time solve for connecting two points, then the optimal unconstrained burn plan. If JuMP won't cooperate, need to write the shooting methods