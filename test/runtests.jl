using DiffEqJumpExtensions
using Base.Test

# write your own tests here
tic()
@time @testset "Linear reaction tests" begin include("linearreaction_test.jl") end
toc()
