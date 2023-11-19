#=
1) Init tree with root = q_init
while q_goal not found:
2) generate random sample
    a) if rand() < p, q_rand = q_goal
    b) else, collisionfreesample
3) find closest node to q_rand, designate q_near
4) generate path from q_near to q_rand
5) take step along that path according to step size r
    a) if collision free, add as q_new
    b) continue stepping until q_rand is reached or collision is found
6) add all collision free paths to tree
7) If (q_new-q_goal).norm() < epsilon, consider goal reached and return
=#

#cpp version has a "findclosestnode" which should already be DONE
#and a "plan" function that did most of the work
#may need some kind of solution struct.