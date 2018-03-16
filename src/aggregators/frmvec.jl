
using FunctionWrappers

mutable struct FRMVECJumpAggregation{T,F1,F2} <: AbstractJumpAggregator
    next_jump_time::T    
    end_time::T
    next_jump::Int64
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
end
  
@inline function (p::FRMVECJumpAggregation)(u, t, integrator) # condition
    p.next_jump_time == t
end
  
@inline function (p::FRMVECJumpAggregation)(integrator) # affect!
    
    @inbounds p.affects![p.next_jump](integrator)

    # find time till next rx, and which rx occurs
    dt,rx = get_next_rx_vec(integrator.u, integrator.p, integrator.t, p.rates)

    # update state
    p.next_jump_time = integrator.t + dt
    p.next_jump      = rx
    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end
    nothing
end

@inline function (p::FRMVECJumpAggregation)(dj, u, t, integrator) # initialize
    # find time till next rx, and which rx occurs
    dt,rx = get_next_rx_vec(u, integrator.p, t, p.rates)

    # update state
    p.next_jump_time = t + dt
    p.next_jump      = rx

    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end

    nothing    
end
  
@inline function get_next_rx_vec(u, params, t, rates)
    dt = typemax(t)
    rx = 0

    for i in eachindex(rates)
        @inbounds dtrx = randexp() / rates[i](u, params, t)
        if dt > dtrx
            dt = dtrx
            rx = i
        end                 
    end

    dt,rx
end

@inline function aggregate(aggregator::FRMVEC, u, p, t, end_time, constant_jumps, save_positions)    
    RateWrapper = FunctionWrappers.FunctionWrapper{Float64,Tuple{typeof(u), typeof(p), typeof(t)}}
    rates = [RateWrapper(c.rate) for c in constant_jumps]
    affects! = [c.affect! for c in constant_jumps]
    dt, rx = get_next_rx_vec(u, p, t, rates)
    
    FRMVECJumpAggregation(t+dt, end_time, rx, rates, affects!, save_positions)
end
  
  DiscreteCallback(c::FRMVECJumpAggregation) = DiscreteCallback(c, c, initialize = c, save_positions = c.save_positions)