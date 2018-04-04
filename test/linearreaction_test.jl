# calculates the mean from N stochastic A->B reactions at different rates
using DiffEqBase
import DiffEqJump, DiffEqJumpExtensions
using Base.Test

doprint  = true
Nsims    = 32000
tf       = .1
baserate = .1
A0       = 100
exactmean = (t,ratevec) -> A0 * exp(-sum(ratevec) * t)


function A_to_B_mean(N, method)
    rates = ones(Float64, N) * baserate;
    cumsum!(rates, rates)    

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
    jset      = DiffEqJumpExtensions.JumpSet((), jumps, nothing, nothing)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, method, jset; save_positions=(false,false))

    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end

    if doprint
        println("samp mean: ", mean(Asamp), ", act mean = ", exactmean(tf, rates))
    end

    @test abs(mean(Asamp) - exactmean(tf, rates)) < 1.
end

function A_to_B_mean_ma(N, method)
    rates = ones(Float64, N) * baserate;
    cumsum!(rates, rates)

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

    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end

    if doprint
        println("samp mean: ", mean(Asamp), ", act mean = ", exactmean(tf, rates))
    end

    @test abs(mean(Asamp) - exactmean(tf, rates)) < 1.
end

function A_to_B_mean_hybrid(N, method)
    rates = ones(Float64, N) * baserate;
    cumsum!(rates, rates)

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
     jumpvec = []
     for i in (switchidx+1):N
         ratefunc = (u,p,t) -> rates[i] * u[1]
         affect!  = function (integrator)
             integrator.u[1] -= 1
             integrator.u[2] += 1
         end
         push!(jumpvec, DiffEqJumpExtensions.ConstantRateJump(ratefunc, affect!))
     end

    jumps     = ((jump for jump in jumpvec)...)
    majumps   = DiffEqJumpExtensions.MassActionJump(rates[1:switchidx] , reactstoch, netstoch)
    jset      = DiffEqJumpExtensions.JumpSet((), jumps, nothing, majumps)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, method, jset; save_positions=(false,false))

    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end

    if doprint
        println("samp mean: ", mean(Asamp), ", act mean = ", exactmean(tf, rates))
    end

    @test abs(mean(Asamp) - exactmean(tf, rates)) < 1.
end



# function wrappers
method = DiffEqJumpExtensions.DirectFunWrappers()
A_to_B_mean(16, method)

# mass action
method = DiffEqJumpExtensions.DirectMassAction()
A_to_B_mean_ma(16, method)

# hybrid
A_to_B_mean_hybrid(16, method)