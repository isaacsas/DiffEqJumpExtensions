using DiffEqBase, OrdinaryDiffEq, DiffEqJump, DiffEqBiological
using BenchmarkTools
using Plots; 

doPlot      = false
doBenchMark = true 

# DNA repression model DiffEqBiological
rs = @reaction_network dtype begin
    k1, DNA --> mRNA + DNA
    k2, mRNA --> mRNA + P
    k3, mRNA --> 0
    k4, P --> 0
    k5, DNA + P --> DNAR
    k6, DNAR --> DNA + P
end k1 k2 k3 k4 k5 k6
tf     = 20000.0
params = (.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.)
prob   = DiscreteProblem([1,0,0,0], (0.0, tf), params)

# SSAs to test
methods = (Direct(),)

# plotting one full trajectory
if doPlot
    plothand = plot()
    for method in methods
        jump_prob = JumpProblem(prob, method, rs)
        sol = solve(jump_prob, SSAStepper())
        plot!(plothand, sol.t, sol[3,:])
    end
    display(plothand)
end

# benchmark performance
if doBenchMark
    # number of save times to look at
    num_save_vec = 25000 #10.^(3:6)

    # exact methods
    for method in methods
        println("Solving with method: ", typeof(method), ", using SSAStepper")
        jump_prob = JumpProblem(prob, method, rs)
        @btime solve($jump_prob, SSAStepper())
    end
    println()

    # methods that save at fixed times
    for num_save in num_save_vec
        dtsave = tf / num_save

        for method in methods    
            println("Solving with method: ", typeof(method), ", using SSAStepper and ", num_save, " points")
            jump_prob = JumpProblem(prob, method, rs, save_positions=(false,false))
            @btime solve($jump_prob, SSAStepper(), saveat=$dtsave)                        
        end 
        
        println()
    end

end

# write your own tests here
#@test 1 == 1
