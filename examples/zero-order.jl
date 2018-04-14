using DiffEqBase, DiffEqJump, Plots
rate = [2.0]
rs = [Vector{Pair{Int,Int}}()]
ns = [[1 => 1]]
jump = MassActionJump(rate, rs, ns)
prob = DiscreteProblem([100],(0.,100.))
jump_prob = JumpProblem(prob, Direct(), jump)
sol = solve(jump_prob, SSAStepper())
plot(sol)
