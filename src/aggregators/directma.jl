
mutable struct DirectMAJumpAggregation{T,S,F1,F2} <: AbstractJumpAggregator
    next_jump_time::T
    next_jump::Int64
    end_time::T
    cur_rates::Vector{T}
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
end

@inline function (p::DirectMAJumpAggregation)(u,t,integrator) # condition
    p.next_jump_time == t
end

function (p::DirectMAJumpAggregation)(integrator) # affect!
    @inbounds executerx!(integrator.u, p.ma_jumps.net_stoch[p.next_jump])
    update_samples(p, integrator, integrator.u, integrator.p, integrator.t)
    nothing
end

function (p::DirectMAJumpAggregation)(dj,u,t,integrator) # initialize
    update_samples(p, integrator, u, integrator.p, t)
    nothing
end

################ DirectMassAction specific implementation functions
function update_samples(p::DirectMAJumpAggregation, integrator, u, params, t)
    
    # next jump time
    p.sum_rate, ttnj = time_to_next_jump_ma(u, params, t, p.ma_jumps, p.cur_rates)
    
    # add jump event to event queue
    @fastmath p.next_jump_time = t + ttnj    
    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end

    # next jump 
    rn = rand() * p.sum_rate
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rn) 
end

@fastmath function time_to_next_jump_ma(u, params, t, majumps::MassActionJump, cur_rates)
    prev_rate = zero(t)
    new_rate  = zero(t)
    @inbounds for i in 1:length(majumps.scaled_rates)        
        new_rate     = evalrxrate(u, majumps.scaled_rates[i], majumps.reactant_stoch[i])
        cur_rates[i] = new_rate + prev_rate
        prev_rate    = cur_rates[i]
    end

    @inbounds sum_rate = cur_rates[end]
    sum_rate, randexp()/sum_rate
end

@inline function aggregate(aggregator::DirectMassAction, u, p, t, end_time,
                            constant_jumps, ma_jumps, save_positions)

    # handle constant jumps using function wrappers
    if constant_jumps != nothing
        rates, affects! = get_jump_info_fwrappers(u,p,t,constant_jumps)
    else
        rates    = []
        affects! = []
    end
    
    # current jump rates, just mass action jumps for now
    cur_rates = Vector{typeof(t)}(length(ma_jumps.scaled_rates))

    DirectMAJumpAggregation(zero(t), 0, end_time, cur_rates, zero(t), 
                            ma_jumps, rates, affects!, save_positions)
end

DiscreteCallback(c::DirectMAJumpAggregation) = DiscreteCallback(c, c, initialize=c, save_positions=c.save_positions)
