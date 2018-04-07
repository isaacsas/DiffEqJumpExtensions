using DiffEqBase, DiffEqJump

using Plots; 
doPlot = true
using BenchmarkTools
doBenchMark = false

Nsims = 8000

# average number of proteins in a simulation
function runSSAs(jump_prob)
    Psamp = zeros(Int, Nsims)
    for i in 1:Nsims
        sol = solve(jump_prob, SSAStepper())
        Psamp[i] = sol[3,end]
    end
    mean(Psamp)
end


# DNA repression model DiffEqBiological
# using DiffEqBiological
# rs = @reaction_network dtype begin
#     k1, DNA --> mRNA + DNA
#     k2, mRNA --> mRNA + P
#     k3, mRNA --> 0
#     k4, P --> 0
#     k5, DNA + P --> DNAR
#     k6, DNAR --> DNA + P
# end k1 k2 k3 k4 k5 k6


# model using mass action jumps
# ids: DNA=1, mRNA = 2, P = 3, DNAR = 4
reactstoch = 
[ 
    [1 => 1],
    [2 => 1],
    [2 => 1],
    [3 => 1],
    [1 => 1, 3 => 1],
    [4 => 1] 
]
netstoch = 
[ 
    [2 => 1],
    [3 => 1],
    [2 => -1],
    [3 => -1],
    [1 => -1, 3 => -1, 4 => 1],
    [1 => 1, 3 => 1, 4 => -1] 
]
rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
majumps = MassActionJumps(rates, reactstoch, netstoch)

tf     = 20000.0
prob   = DiscreteProblem([1,0,0,0], (0.0, tf), rates)

# SSAs to test
SSAalgs = (Direct(), DirectFW(), FRM(), FRMFW())

# plotting one full trajectory
if doPlot
    plothand = plot()
    for alg in SSAalgs
        jump_prob = JumpProblem(prob, alg, majumps)
        sol = solve(jump_prob, SSAStepper())
        plot!(plothand, sol.t, sol[3,:])
    end
    display(plothand)
end

# benchmark performance
if doBenchMark
    # exact methods
    for alg in SSAalgs
        println("Solving with method: ", typeof(method), ", using SSAStepper")
        jump_prob = JumpProblem(prob, method, rs)
        @btime solve($jump_prob, SSAStepper())
    end
    println()
end

# write your own tests here
#@test 1 == 1
