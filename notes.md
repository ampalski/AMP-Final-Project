General flow:
* Sample n points
* do pre-computation on unbounded connections with limiting delta-v and tof
* Run through the base technique, save off results
* Send results through smoothing process, save off results
* Compare results... accuracy, delta-v, runtimes, etc.

Data structures:
* Graph holds which nodes are connected (directed, so children are connections)
* liveGraph holds connections as the planner makes them, directed
* edgeCosts sparse array holds delta-v cost for edge costs
* fullEdgeCosts dict (int, int) -> (vector, vector, float) to store dv1, dv2, T for connections ... 
* Map holds int -> vector{int} for which other nodes are in its neighborhood after pre-processing -- hold off, I think the graph itself can do that
* samples Dict holds sample locations (int, node number -> vector RIC/RICdot)
* Simple array (vector) holds the Vertex_open, Vertex_unvisited, other sets. Just need to store node number as int -- implemented as sets for now for fast checking/insertion/deletion
* Probably have a struct hold all the common features, separate structs for any technique specific data
* Primary waypoint paths can just come back as matrices. But want the ability to get matrices for every path held in the graph for plotting

Things I need to figure out:
* Linear programming. Make sure I can do the min-time solve for connecting two points, then the optimal unconstrained burn plan. If JuMP won't cooperate, need to write the shooting methods

Ideas for paper:
* Plot of # of nodes vs success rate for RRT. It seems like numNodes is currently a limiting factor, especially with goal bias set to 0.05
* 2000 samples took 17746s to connect 26992.8
