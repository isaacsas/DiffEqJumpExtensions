module DiffEqJumpExtensions

using DiffEqBase, Compat, Requires
using DiffEqJump

import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices
import DiffEqJump: aggregate, get_jump_info_tuples, get_jump_info_fwrappers, executerx!, 
                   evalrxrate, fill_cur_rates, register_next_jump_time!, build_jump_aggregation
import Base: size, getindex, setindex!, length, similar, indices, show

#include("massaction_rates.jl")
include("aggregators/aggregators.jl")
include("aggregators/frm.jl")

export FRM, FRMFW, FRMJumpAggregation, aggregate

end # module
