using DiffEqBase, DiffEqJump
using Plots; plotlyjs()

tf      = .01
u0      = [200, 100, 150]

# MODEL SETUP

# DiffEqBiological Version:
using DiffEqBiological
rs = @reaction_network dtype begin
    k1, 2A --> B
    k2, B --> 2A 
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
rates = [1., 2., .5, .75, .25]

# model using mass action jumps
# ids: A = 1, B = 2, C = 3

# reactant stoichiometry is stored as a Vector{Vector{Pairs}}
# the kth entry in reactstoch is the stoichiometry of the kth reaction
reactstoch = 
[ 
    [1 => 2],            # 2A
    [2 => 1],            # B
    [1 => 1, 2 => 1],    # A, B
    [3 => 1],            # C
    [3 => 3]             # 3C
]

# net stoichiometry is stored as a Vector{Vector{Pairs}}
# the kth entry in netstoch is the net stoichiometry of the kth reaction
netstoch = 
[ 
    [1 => -2, 2 => 1],          # A -= 2, B += 1
    [1 => 2, 2 => -1],          # A += 2, B -= 1
    [1 => -1, 2 => -1, 3 => 1], # A -= 1, B -= 1, C += 1
    [1 => 1, 2 => 1, 3 => -1],  # A += 1, B += 1, C -= 1
    [1 => 3, 3 => -3]           # A += 3, C -= 3
]

# how to create the mass action jump type
majumps = MassActionJump(rates, reactstoch, netstoch)


prob = DiscreteProblem(u0, (0.0, tf), rates)

# solving with DiffEqBiological
jump_prob = JumpProblem(prob, Direct(), rs)
sol = solve(jump_prob, SSAStepper())
plothand = plot(sol, seriestype=:steppost, reuse=false)
display(plothand)


# solving with majumps
jump_prob = JumpProblem(prob, Direct(), majumps)
sol = solve(jump_prob, SSAStepper())
plothand = plot(sol, seriestype=:steppost, reuse=false)
display(plothand)

