using DiffEqBase, DiffEqJump, DiffEqBiological
using DiffEqJumpExtensions
using BenchmarkTools
using BioSimulator


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
tf = 200.
params = (1.,)
prob   = DiscreteProblem(10*ones(Int64,N), (0.0, tf), params)

# Biosimlator model of the same network
model = Network("Diffusion")
for i in 1:N
    model <= Species("X$(i)", 10)
end
for i in 1:(N-1)
    model <= Reaction("diffusion$(i)a", params[1], "X$(i) --> X$(i+1)")
    model <= Reaction("diffusion$(i)b", params[1], "X$(i+1) --> X$(i)")
end

# mass action version
scaled_rates   = ones(Float64,2*(N-1))
reactant_stoch = Vector{Vector{Pair{Int64,Int64}}}();
net_stoch      = Vector{Vector{Pair{Int64,Int64}}}();
for i = 1:(N-1)
    push!(reactant_stoch, [i => 1])
    push!(reactant_stoch, [(i+1) => 1])
    push!(net_stoch, [i => -1, (i+1) => 1])
    push!(net_stoch, [(i+1) => -1, i => 1])
end
# reactant_stoch = Vector{Array{Int64,2}}();
# net_stoch      = Vector{Array{Int64,2}}();
# for i = 1:(N-1)
#     push!(reactant_stoch, [i,1]')
#     push!(reactant_stoch, [(i+1),1]')
#     push!(net_stoch, [i (i+1); -1 1])
#     push!(net_stoch, [(i+1) i; -1  1])
# end

# DiffEqJump methods to benchmark
#methods = (Direct2(), DirectVEC(), DiffEqJumpExtensions.FRM(), FRMVEC())
#methods = (DirectVEC(), FRMVEC())
methods = (DirectMA(scaled_rates, reactant_stoch, net_stoch), DirectVEC())


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
