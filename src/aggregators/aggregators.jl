
struct DirectFunWrappers <: AbstractAggregatorAlgorithm end
struct DirectMassAction <: AbstractAggregatorAlgorithm end

abstract type AbstractSSAJumpAggregator <: AbstractJumpAggregator end 

DiscreteCallback(c::AbstractSSAJumpAggregator) = DiscreteCallback(c, c, initialize=c, save_positions=c.save_positions)

# tells integrator when the next jump will occur
function register_next_jump_time!(integrator, p::AbstractSSAJumpAggregator, t)        
    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end
end

