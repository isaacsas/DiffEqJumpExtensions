module DiffEqJumpExtensions

using DiffEqBase, Compat, Requires
using DiffEqJump

import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices
import Base: size, getindex, setindex!, length, similar, indices, show

include("massaction_rates.jl")
include("jumps.jl")
include("aggregators/aggregators.jl")
include("aggregators/direct.jl")
include("aggregators/directma.jl")


export JumpProblem

export DirectMA


end # module
