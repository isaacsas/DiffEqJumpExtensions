using DiffEqBase, DiffEqJump, DiffEqBiological
using BenchmarkTools

doBenchmark = true

# number of save times to look at
num_save_vec = 25000 

# diffusion model
N = 32
println("N = ", N)
diffnetwork = "@reaction_network difftype begin\n"
for i in 1:(N-1)
    diffnetwork *= "\t K, X$(i) --> X$(i+1)\n"
    diffnetwork *= "\t K, X$(i+1) --> X$(i)\n"
end
diffnetwork *= "end K"
rs = eval( parse(diffnetwork) )
tf = 10.
params = (1.,)
prob   = DiscreteProblem(10*ones(Int64,N), (0.0, tf), params)

methods = (Direct(),)

if doBenchmark
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

        println()
    end
end