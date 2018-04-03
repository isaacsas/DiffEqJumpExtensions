# N stochastic A->B reactions at different rates
using DiffEqBase, DiffEqJump, DiffEqJumpExtensions
using Base.Test

doprint = true

function A_to_B_mean(N)
    Nsims    = 32000
    tf       = .1
    baserate = .1
    A0       = 100

    rates    = ones(Float64, N) * baserate;
    cumsum!(rates, rates)
    exactmean = (t) -> A0*exp(-sum(rates) * t)

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
    jumps = ((jump for jump in jumpvec)...)
    jset  = DiffEqJumpExtensions.JumpSet((), jumps, nothing, nothing)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, DiffEqJumpExtensions.DirectFunWrappers(), jset; save_positions=(false,false))

    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end

    if doprint
        println("samp mean: ", mean(Asamp), ", act mean = ", exactmean(tf))
    end

    @test abs(mean(Asamp) - exactmean(tf)) < 1.

end

function A_to_B_mean_ma(N)
    Nsims    = 32000
    tf       = .1
    baserate = .1
    A0       = 100

    rates = ones(Float64, N) * baserate;
    cumsum!(rates, rates)
    exactmean = (t) -> A0*exp(-sum(rates) * t)

    reactstoch = Vector{Vector{Pair{Int64,Int64}}}();
    netstoch   = Vector{Vector{Pair{Int64,Int64}}}();
    for i = 1:N
        push!(reactstoch,[1 => 1])
        push!(netstoch,[1 => -1, 2=>1])
    end


    majumps   = DiffEqJumpExtensions.MassActionJump(rates, reactstoch, netstoch)
    jset      = DiffEqJumpExtensions.JumpSet((), (), nothing, majumps)
    prob      = DiscreteProblem([A0,0], (0.0,tf))
    jump_prob = DiffEqJumpExtensions.JumpProblem(prob, DiffEqJumpExtensions.DirectMassAction(), jset; save_positions=(false,false))

    Asamp = zeros(Int64,Nsims)
    for i in 1:Nsims
        sol = DiffEqJumpExtensions.solve(jump_prob, DiffEqJumpExtensions.SSAStepper())
        Asamp[i] = sol[1,end]
    end

    if doprint
        println("samp mean: ", mean(Asamp), ", act mean = ", exactmean(tf))
    end

    @test abs(mean(Asamp) - exactmean(tf)) < 1.

end



# tuples
A_to_B_mean(5)

# function wrappers
A_to_B_mean(15)

# mass action
A_to_B_mean_ma(15)