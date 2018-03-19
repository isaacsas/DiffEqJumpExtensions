using FunctionWrappers

mutable struct DirectMAJumpAggregation{T,S1,S2,F1,F2} <: AbstractJumpAggregator
    next_jump::T
    end_time::T
    cur_rates::Vector{T}
    sum_rate::T
    aggregator::DirectMA{S1,S2}
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
end

@inline function (p::DirectMAJumpAggregation)(u,t,integrator) # condition
    p.next_jump==t
end

function (p::DirectMAJumpAggregation)(integrator) # affect!
    rng_val = rand() * p.sum_rate
    i       = searchsortedfirst(p.cur_rates, rng_val)
    @inbounds executerx!(integrator.u, p.aggregator.net_stoch[i])
    p.sum_rate,ttnj = time_to_next_jump_ma(integrator.u,integrator.p,integrator.t,p.aggregator,p.cur_rates)
    @fastmath p.next_jump = integrator.t + ttnj
    if p.next_jump < p.end_time
        add_tstop!(integrator,p.next_jump)
    end
    nothing
end

function (p::DirectMAJumpAggregation)(dj,u,t,integrator) # initialize
    sum_rate,next_jump = time_to_next_jump_ma(u,integrator.p,t,p.aggregator,p.cur_rates)
    p.sum_rate         = sum_rate
    p.next_jump        = t + next_jump
    if p.next_jump < p.end_time
        add_tstop!(integrator,p.next_jump)
    end
    nothing
end

@fastmath function time_to_next_jump_ma(u,p,t,aggregator::DirectMA,cur_rates)
    @inbounds sum_rate     = evalrxrate(u, aggregator.scaled_rates[1], aggregator.reactant_stoch[1])
    @inbounds cur_rates[1] = sum_rate
    @inbounds for i in 2:length(cur_rates)
        new_rate     = evalrxrate(u, aggregator.scaled_rates[i], aggregator.reactant_stoch[i])
        sum_rate    += new_rate
        cur_rates[i] = new_rate + cur_rates[i-1]
    end

    sum_rate,randexp()/sum_rate
end

@inline function aggregate(aggregator::DirectMA,u,p,t,end_time,constant_jumps,save_positions)
    # RateWrapper        = FunctionWrappers.FunctionWrapper{typeof(t),Tuple{typeof(u), typeof(p), typeof(t)}}
    # rates              = [RateWrapper(c.rate) for c in constant_jumps]
    # AffectWrapper      = FunctionWrappers.FunctionWrapper{Void,Tuple{Any}}
    # affects!           = [AffectWrapper(x->(c.affect!(x);nothing)) for c in constant_jumps]
    rates = [evalrxrate]
    affects! = [executerx!]
    cur_rates          = Vector{typeof(t)}(length(aggregator.scaled_rates))
    sum_rate,next_jump = time_to_next_jump_ma(u,p,t,aggregator,cur_rates)
    DirectMAJumpAggregation(next_jump, end_time, cur_rates, sum_rate, aggregator,
                            rates, affects!, save_positions)
end

DiscreteCallback(c::DirectMAJumpAggregation) = DiscreteCallback(c,c,initialize=c,save_positions=c.save_positions)
