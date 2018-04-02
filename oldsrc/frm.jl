
mutable struct FRMJumpAggregation{T,F1,F2} <: AbstractJumpAggregator
    next_jump_time::T    
    end_time::T
    next_jump::Int64
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
end
  
@inline function (p::FRMJumpAggregation)(u, t, integrator) # condition
    p.next_jump_time == t
end
  
@inline function (p::FRMJumpAggregation)(integrator) # affect!
    
    @inbounds p.affects![p.next_jump](integrator)

    # find time till next rx, and which rx occurs
    dt,rx = get_next_rx(integrator.u, integrator.p, integrator.t, p.rates)

    # update state
    p.next_jump_time = integrator.t + dt
    p.next_jump      = rx
    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end
    nothing
end

@inline function (p::FRMJumpAggregation)(dj, u, t, integrator) # initialize
    # find time till next rx, and which rx occurs
    dt,rx = get_next_rx(u, integrator.p, t, p.rates)

    # update state
    p.next_jump_time = t + dt
    p.next_jump      = rx

    if p.next_jump_time < p.end_time
        add_tstop!(integrator, p.next_jump_time)
    end

    nothing    
end
  
@inline function get_next_rx(u, params, t, rates)
    dt = typemax(t)
    rx = 0

    dt,rx = get_next_rx(u, params, t, 1, rates...)

    #dt,rx
end

@inline function get_next_rx(u, params, t, idx, rate, rates...)
    dt2   = randexp() / rate(u, params, t)
    dt,rx = get_next_rx(u, params, t, idx+1, rates...)

    if dt2 < dt
        dt = dt2
        rx = idx
    end
        
    dt,rx
end

@inline function get_next_rx(u, params, t, idx, rate)
    dt = randexp() / rate(u, params, t)
    dt,idx
end

@inline function aggregate(aggregator::FRM, u, p, t, end_time, constant_jumps, save_positions)
    rates = ((c.rate for c in constant_jumps)...)
    affects! = ((c.affect! for c in constant_jumps)...)
    dt, rx = get_next_rx(u, p, t, rates)
    
    FRMJumpAggregation(t+dt, end_time, rx, rates, affects!, save_positions)
end
  
  DiscreteCallback(c::FRMJumpAggregation) = DiscreteCallback(c, c, initialize = c, save_positions = c.save_positions)