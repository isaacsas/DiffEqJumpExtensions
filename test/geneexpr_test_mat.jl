using DiffEqBase, DiffEqJump, DiffEqJumpExtensions
using Base.Test

using Plots; plotlyjs()
doplot = true
using BenchmarkTools
dobenchmark = true

dotestmean   = false
doprintmeans = false

# SSAs to test
#SSAalgs = (SortingDirect(),Direct()) #, DirectFW(), FRM(), FRMFW())
SSAalgs = (Direct(),)

# Nsims        = 80000
# tf           = 100000.0
Nsims        = 100
tf           = 10000.0
u0           = [1,0,0,0]
expected_avg = 5.926553750000000e+02
reltol       = .01

# average number of proteins in a simulation
function runSSAs(jump_prob)
    Psamp = zeros(Int, Nsims)
    for i in 1:Nsims
        sol = solve(jump_prob, SSAStepper())
        Psamp[i] = sol[3,end]
    end
    mean(Psamp)
end

# MODEL SETUP

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
    [1 1]',
    [2 1]',
    [2 1]',
    [3 1]',
    [1 1; 3 1]',
    [4 1]'
]
netstoch =
[
    [2 1]',
    [3 1]',
    [2 -1]',
    [3 -1]',
    [1 -1; 3 -1; 4 1]',
    [1 1; 3 1; 4 -1]'
]
rates = [50, 50*(20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
#rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
majumps = MassActionJump(rates, reactstoch, netstoch)


# TESTING:
prob = DiscreteProblem(u0, (0.0, tf), rates)

# plotting one full trajectory
if doplot
    plothand = plot(reuse=false)
    for alg in SSAalgs
        jump_prob = JumpProblem(prob, alg, majumps, save_positions=(false,false))
        sol = solve(jump_prob, SSAStepper(), saveat=(tf/1000.))
        plot!(plothand, sol.t, sol[3,:], seriestype=:steppost)
    end
    display(plothand)
end

# test the means
if dotestmean
    means = zeros(Float64,length(SSAalgs))
    for (i,alg) in enumerate(SSAalgs)
        jump_prob = JumpProblem(prob, alg, majumps, save_positions=(false,false))
        means[i]  = runSSAs(jump_prob)
        relerr = abs(means[i] - expected_avg) / expected_avg
        if doprintmeans
            println("Mean from method: ", typeof(alg), " is = ", means[i], ", rel err = ", relerr)
        end

        if dobenchmark
            @btime (runSSAs($jump_prob);)
        end


        @test abs(means[i] - expected_avg) < reltol*expected_avg
    end
end


# benchmark performance
if dobenchmark
    # exact methods
    for alg in SSAalgs
        println("Solving with method: ", typeof(alg), ", using SSAStepper")
        jump_prob = JumpProblem(prob, alg, majumps, save_positions=(false,false))
        @btime solve($jump_prob, SSAStepper(), saveat=(tf/1000.))
        @btime (runSSAs($jump_prob);)
    end
    println()
end
