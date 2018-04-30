module DiffEqJumpExtensions

using DiffEqBase, Compat, Requires, DataStructures
using DiffEqJump

import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices
import DiffEqJump: aggregate, get_jump_info_tuples, get_jump_info_fwrappers, executerx!, 
                   evalrxrate, fill_cur_rates, register_next_jump_time!, build_jump_aggregation, 
                   get_num_majumps, make_dependency_graph
import Base: size, getindex, setindex!, length, similar, indices, show

include("aggregators/aggregators.jl")
include("aggregators/nrm.jl")

export NRM, aggregate

end # module
