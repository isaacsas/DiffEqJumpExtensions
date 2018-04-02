using DiffEqJumpExtensions
using Base.Test

# write your own tests here
tic()
@time @testset "Functionwrapper tests" begin include("functwrapper_tests.jl") end
toc()
