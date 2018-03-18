using DiffEqBase, DiffEqJump, DiffEqBiological
using DiffEqJumpExtensions
using BenchmarkTools
using BioSimulator

# DiffEqJump methods to benchmark
#methods = (Direct2(), DirectVEC(), DiffEqJumpExtensions.FRM(), FRMVEC())
#methods = (DirectVEC(), FRMVEC())
methods = (DirectVEC(),)

# number of save times to look at
num_save_vec = 25000 #10.^(3:6)

# diffusion model
N = 256
diffnetwork = "@reaction_network difftype begin\n"
for i in 1:(N-1)
    diffnetwork *= "\t K, X$(i) --> X$(i+1)\n"
    diffnetwork *= "\t K, X$(i+1) --> X$(i)\n"
end
diffnetwork *= "end K"
#println(diffnetwork)
rs = eval( parse(diffnetwork) )
tf = 20.
params = (1.,)
prob   = DiscreteProblem(10*ones(N), (0.0, tf), params)

# Biosimlator model of the same network
model = Network("Diffusion")
for i in 1:N
    model <= Species("X$(i)")
end
for i in 1:(N-1)
    model <= Reaction("diffusion$(i)a", params[1], "X$(i) --> X$(i+1)")
    model <= Reaction("diffusion$(i)b", params[1], "X$(i+1) --> X$(i)")
end


# exact methods
# println("Full solution paths saving:\n")
# for method in methods
#     println("Solving with method: ", typeof(method), ", using SSAStepper")
#     jump_prob = JumpProblem(prob, method, rs)
#     @btime solve($jump_prob, SSAStepper())
# end
# println()

# methods with different numbers of points
println("Saving at fixed times:\n")
for num_save in num_save_vec
    dtsave = tf / num_save

    for method in methods
        println("Solving with method: ", typeof(method), ", using SSAStepper and ", num_save, " points")
        jump_prob = JumpProblem(prob, method, rs, save_positions=(false,false))
        sol = solve(jump_prob, SSAStepper(), saveat=dtsave);
        @btime solve($jump_prob, SSAStepper(), saveat=$dtsave)
    end

    # BioSimulator:
    println("BioSimulator, Direct, saving at ", num_save, " points:")
    result = simulate(model, BioSimulator.SSA, time=tf, epochs=num_save, trials=1);
    @btime simulate($model, BioSimulator.SSA, time=$tf, epochs=$num_save, trials=1);

    # println("BioSimulator, FRM, saving at ", num_save, " points:")
    # @btime simulate($model, BioSimulator.FRM, time=$tf, epochs=$num_save, trials=1);

    println()
end
