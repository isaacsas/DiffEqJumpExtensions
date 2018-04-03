
struct DirectFunWrappers <: AbstractAggregatorAlgorithm end

# SSAs subtype the AbstractSSAAggregatorAlgorithm
# all AbstractSSAAggregators are assumed to have the following fields:
# next_jump      = index of the next jump
# next_jump_time = time of the next jump
# end_time       = time to stop the simulation
# rng            = random number generator
# save_positions::Tuple{Bool,Bool} = when to save the simulation data at
abstract type AbstractSSAAggregatorAlgorithm <: AbstractAggregatorAlgorithm end
abstract type AbstractSSAJumpAggregator <: AbstractJumpAggregator end 
DiscreteCallback(c::AbstractSSAJumpAggregator) = DiscreteCallback(c, c, initialize=c, save_positions=c.save_positions)

# tells integrator when the next jump will occur
function register_next_jump_time!(integrator, p::AbstractSSAJumpAggregator, t)        
    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end
end

struct DirectMassAction <: AbstractSSAAggregatorAlgorithm end



