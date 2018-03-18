module DiffEqJumpExtensions

using DiffEqBase, Compat, Requires
using DiffEqJump, DiffEqBiological

import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices
import Base: size, getindex, setindex!, length, similar, indices, show

include("aggregators/aggregators.jl")
include("aggregators/direct2.jl")
include("aggregators/directvec.jl")
include("aggregators/directma.jl")
include("aggregators/frm.jl")
include("aggregators/frmvec.jl")
include("problem.jl")
include("SSA_stepper.jl")
include("massaction_rates.jl")

export JumpProblem

export FRM, FRMVEC, Direct2, DirectVEC, DirectMA


end # module
