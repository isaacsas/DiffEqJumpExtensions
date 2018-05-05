# from McCollum et al, 
# "The sorting direct method for stochastic simulation of biochemical systems with varying reaction execution behavior"
# Comp. Bio. and Chem., 30, pg. 39-49 (2006). 

using DiffEqBase, DiffEqJump, DiffEqBiological
using Base.Test

using Plots; plotlyjs()
doplot = false
using BenchmarkTools
dobenchmark = true

dotestmean   = false
doprintmeans = false

# SSAs to test
#SSAalgs = (SortingDirect(),Direct()) #, DirectFW(), FRM(), FRMFW())
SSAalgs = (Direct(), SortingDirect(),)

# Nsims        = 80000
# tf           = 100000.0
Nsims        = 100
tf           = 500.0
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
N = 10
genenetwork = "@reaction_network gtype begin\n"
for i in 1:N
    genenetwork *= "\t 10.0, G$(2*i-1) --> G$(2*i-1) + M$(2*i-1)\n"
    genenetwork *= "\t 10.0, M$(2*i-1) --> M$(2*i-1) + P$(2*i-1)\n"
    genenetwork *= "\t 1.0,  M$(2*i-1) --> 0\n"
    genenetwork *= "\t 1.0,  P$(2*i-1) --> 0\n"

    genenetwork *= "\t 5.0, G$(2*i) --> G$(2*i) + M$(2*i)\n"
    genenetwork *= "\t 5.0, M$(2*i) --> M$(2*i) + P$(2*i)\n"
    genenetwork *= "\t 1.0,  M$(2*i) --> 0\n"
    genenetwork *= "\t 1.0,  P$(2*i) --> 0\n"

    genenetwork *= "\t 0.0001, G$(2*i) + P$(2*i-1) --> G$(2*i)_ind \n"
    genenetwork *= "\t 100., G$(2*i)_ind --> G$(2*i)_ind + M$(2*i)\n"
end
genenetwork *= "end"
rs = eval( parse(genenetwork) )
u0 = zeros(Int, length(rs.syms))
for i = 1:(2*N)
    u0[findfirst(rs.syms, Symbol("G$(i)"))] = 1    
end

# TESTING:
prob = DiscreteProblem(u0, (0.0, tf))

# plotting one full trajectory
if doplot
    plothand = plot(reuse=false)
    for alg in SSAalgs
        jump_prob = JumpProblem(prob, alg, rs, save_positions=(false,false))
        sol = solve(jump_prob, SSAStepper(), saveat=(tf/1000.))
        plot!(plothand, sol.t, sol[[3,6],:]', lab=string(alg))
    end
    display(plothand)
end

# # test the means
# if dotestmean
#     means = zeros(Float64,length(SSAalgs))
#     for (i,alg) in enumerate(SSAalgs)
#         jump_prob = JumpProblem(prob, alg, majumps, save_positions=(false,false))
#         means[i]  = runSSAs(jump_prob)
#         relerr = abs(means[i] - expected_avg) / expected_avg
#         if doprintmeans
#             println("Mean from method: ", typeof(alg), " is = ", means[i], ", rel err = ", relerr)
#         end

#         if dobenchmark
#             @btime (runSSAs($jump_prob);)
#         end


#         @test abs(means[i] - expected_avg) < reltol*expected_avg
#     end
# end


# # benchmark performance
if dobenchmark
    # exact methods
    for alg in SSAalgs
        println("Solving with method: ", typeof(alg), ", using SSAStepper")
        jump_prob = JumpProblem(prob, alg, rs, save_positions=(false,false))
        @btime solve($jump_prob, SSAStepper(), saveat=(tf/1000.))
        #@btime (runSSAs($jump_prob);)
    end
    println()
end
