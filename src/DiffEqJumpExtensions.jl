module DiffEqJumpExtensions

using DiffEqBase, Compat, Requires

import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices
import DiffEqJump: AbstractJump, AbstractAggregatorAlgorithm, AbstractJumpAggregator, AbstractJumpProblem
import Base: size, getindex, setindex!, length, similar, indices, show

include("massaction_rates.jl")
include("jumps.jl")
include("problem.jl")
include("aggregators/aggregators.jl")
include("aggregators/directfwrap.jl")
include("aggregators/directma.jl")
include("SSA_stepper.jl")

export MassActionJump

export JumpProblem

export DirectFunWrappers, DirectMassAction


end # module
