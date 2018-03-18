using DiffEqBase, OrdinaryDiffEq, DiffEqJump, DiffEqBiological
using StaticArrays, DataStructures
using DiffEqJumpExtensions
using Plots; pyplot()
using BenchmarkTools
using BioSimulator

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
tf        = 20000.0
params    =  (.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.)
prob      = DiscreteProblem([1,0,0,0], (0.0, tf), params)

# Biosimlator model of the same network
model = Network("Gene Expression with Neg Feedback")
model <= Species("DNA", 1)
model <= Species("mRNA", 0)
model <= Species("P", 0)
model <= Species("DNAR", 0)
model <= Reaction("transcription", .5, "DNA --> mRNA + DNA")
model <= Reaction("translation", 20*log(2.)/120., "mRNA --> mRNA + P")
model <= Reaction("mRNA degradation", (log(2.)/120.), "mRNA --> 0")
model <= Reaction("protein degradation", (log(2.)/600.), "P --> 0")
model <= Reaction("DNA repression", .025, "DNA + P --> DNAR")
model <= Reaction("DNA freedom", 1., "DNAR --> DNA + P")


spec = OrderedDict( :DNA => 1, :mRNA => 2, :P => 3, :DNAR => 4)
scaled_rates = [rate for rate in params]
reactant_stoch = [
    SMatrix{2,1}([spec[:DNA] 1]),
    SMatrix{2,1}([spec[:mRNA] 1]),
    SMatrix{2,1}([spec[:mRNA] 1]),
    SMatrix{2,1}([spec[:P] 1]),
    SMatrix{2,2}([spec[:DNA] spec[:P]; 1 1]),
    SMatrix{2,1}([spec[:DNAR] 1])
    ]
net_stoch = [
       SMatrix{2,1}([spec[:mRNA] 1]'),
       SMatrix{2,1}([spec[:P] 1]'),
       SMatrix{2,1}([spec[:mRNA] -1]'),
       SMatrix{2,1}([spec[:P] -1]'),
       SMatrix{2,3}([spec[:DNA] spec[:P] spec[:DNAR]; -1 -1 1]),
       SMatrix{2,3}([spec[:DNA] spec[:P] spec[:DNAR]; 1 1 -1]) 
       ]


# DiffEqJump methods to benchmark
#methods = (Direct2(), DirectVEC(), DiffEqJumpExtensions.FRM(), FRMVEC())
methods = (DirectVEC(), DirectMA(scaled_rates, reactant_stoch, net_stoch))

# plotting solutions
if doPlot
    result  = simulate(model, algorithm=BioSimulator.SSA, time=tf, epochs=25000, trials=1)
    gui( plot(result.t_index, result[1,:,1], reuse=false) )

    for method in methods
        jump_prob = JumpProblem(prob, method, rs)
        sol = solve(jump_prob, SSAStepper())
        gui(plot!(sol.t, sol[3,:]))
    end

end

if doBenchMark
    # number of save times to look at
    num_save_vec = 25000 #10.^(3:6)

    # exact methods
    for method in methods
        # println("Solving with method: ", typeof(method), ", using FunctionMap")
        # jump_prob = JumpProblem(prob, method, rs)    
        # @btime solve($jump_prob, FunctionMap())

        println("Solving with method: ", typeof(method), ", using SSAStepper")
        jump_prob = JumpProblem(prob, method, rs)
        @btime solve($jump_prob, SSAStepper())
    end
    println()

    # methods with different numbers of points
    for num_save in num_save_vec
        dtsave = tf / num_save

        for method in methods    
            # println("Solving with method: ", typeof(method), ", using FunctionMap and ", num_save, " points")
            # jump_prob = JumpProblem(prob, method, rs, save_positions=(false,false))    
            # @btime solve($jump_prob, FunctionMap(), saveat=$dtsave)                

            println("Solving with method: ", typeof(method), ", using SSAStepper and ", num_save, " points")
            jump_prob = JumpProblem(prob, method, rs, save_positions=(false,false))
            @btime solve($jump_prob, SSAStepper(), saveat=$dtsave)                        
        end 

        # BioSimulator:
        println("BioSimulator, Direct, saving at ", num_save, " points:")
        @btime simulate($model, BioSimulator.SSA, time=$tf, epochs=$num_save, trials=1);
        
        println("BioSimulator, FRM, saving at ", num_save, " points:")
        @btime simulate($model, BioSimulator.FRM, time=$tf, epochs=$num_save, trials=1);
        
        println()
    end

end

# write your own tests here
#@test 1 == 1
