# calculates the mean from N stochastic A->B reactions at different rates
using DiffEqBase
import DiffEqJump, DiffEqJumpExtensions
using Base.Test
using BenchmarkTools

doprint     = true
dobenchmark = true
dotest      = true
Nrxs        = 16
Nsims       = 32000
tf          = .1
baserate    = .1
A0          = 100
exactmean   = (t,ratevec) -> A0 * exp(-sum(ratevec) * t)

rates = ones(Float64, Nrxs) * baserate;
cumsum!(rates, rates)    

function runSSAs(jump_prob)
    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end
    mean(Asamp)
end

function A_to_B_mean_orig(N, method)
    # jump reactions
    jumpvec = []
    for i in 1:N
        ratefunc = (u,p,t) -> rates[i] * u[1]
        affect!  = function (integrator)
            integrator.u[1] -= 1
            integrator.u[2] += 1
        end
        push!(jumpvec, DiffEqJumpExtensions.ConstantRateJump(ratefunc, affect!))
    end

    # convert jumpvec to tuple to send to JumpProblem...
    jumps     = ((jump for jump in jumpvec)...)
    jset      = DiffEqJump.JumpSet((), jumps, nothing)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJump.JumpProblem(prob, method, jset; save_positions=(false,false))

    jump_prob
end

function A_to_B_mean(N, method)
    # jump reactions
    jumps = Vector{DiffEqJumpExtensions.ConstantRateJump}()
    for i in 1:N
        ratefunc = (u,p,t) -> rates[i] * u[1]
        affect!  = function (integrator)
            integrator.u[1] -= 1
            integrator.u[2] += 1
        end
        push!(jumps, DiffEqJumpExtensions.ConstantRateJump(ratefunc, affect!))
    end

    # convert jumpvec to tuple to send to JumpProblem...
    jset      = DiffEqJumpExtensions.JumpSet((), jumps, nothing, nothing)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, method, jset; save_positions=(false,false))

    jump_prob
end

function A_to_B_mean_ma(N, method)
    reactstoch = Vector{Vector{Pair{Int64,Int64}}}();
    netstoch   = Vector{Vector{Pair{Int64,Int64}}}();
    for i = 1:N
        push!(reactstoch,[1 => 1])
        push!(netstoch,[1 => -1, 2=>1])
    end

    majumps   = DiffEqJumpExtensions.MassActionJump(rates, reactstoch, netstoch)
    jset      = DiffEqJumpExtensions.JumpSet((), (), nothing, majumps)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, method, jset; save_positions=(false,false))

    jump_prob
end

function A_to_B_mean_hybrid(N, method)
    # half reactions are treated as mass action and half as constant jumps
    switchidx = (N//2).num

    # mass action reactions
    reactstoch = Vector{Vector{Pair{Int64,Int64}}}();
    netstoch   = Vector{Vector{Pair{Int64,Int64}}}();
    for i in 1:switchidx
        push!(reactstoch,[1 => 1])
        push!(netstoch,[1 => -1, 2=>1])
    end

     # jump reactions
     jumps = Vector{DiffEqJumpExtensions.ConstantRateJump}()
     for i in (switchidx+1):N
         ratefunc = (u,p,t) -> rates[i] * u[1]
         affect!  = function (integrator)
             integrator.u[1] -= 1
             integrator.u[2] += 1
         end
         push!(jumps, DiffEqJumpExtensions.ConstantRateJump(ratefunc, affect!))
     end

    majumps   = DiffEqJumpExtensions.MassActionJump(rates[1:switchidx] , reactstoch, netstoch)
    jset      = DiffEqJumpExtensions.JumpSet((), jumps, nothing, majumps)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, method, jset; save_positions=(false,false))

    jump_prob
end


means = []

# tuples
# method = DiffEqJump.Direct()
# jump_prob_orig = A_to_B_mean_orig(Nrxs, method)
# push!(means, runSSAs(jump_prob_orig))

# function wrappers
method = DiffEqJumpExtensions.DirectMassAction()
jump_prob_fw = A_to_B_mean(Nrxs, method)
push!(means, runSSAs(jump_prob_fw))

# mass action
jump_prob_ma = A_to_B_mean_ma(Nrxs, method)
push!(means, runSSAs(jump_prob_ma))

# hybrid
jump_prob_hybrid = A_to_B_mean_hybrid(Nrxs, method)
push!(means, runSSAs(jump_prob_hybrid))

exactmeanval = exactmean(tf, rates)
for meanval in means
    if doprint
        println("samp mean: ", meanval, ", act mean = ", exactmeanval)
    end

    if dotest
        @test abs(meanval - exactmeanval) < 1.
    end
end

if dobenchmark
    # @btime runSSAs($jump_prob_orig)
    @btime runSSAs($jump_prob_fw)
    @btime runSSAs($jump_prob_ma)
    @btime runSSAs($jump_prob_hybrid)
end



# agg = jump_prob_fw.jump_callback.discrete_callbacks[1].affect!
# integrator = init(jump_prob_fw, DiffEqJumpExtensions.SSAStepper())
# initialize!(agg,integrator,integrator.u, integrator.p,integrator.t)