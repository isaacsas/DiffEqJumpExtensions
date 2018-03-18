using FunctionWrappers

mutable struct DirectVECJumpAggregation{T,F1,F2} <: AbstractJumpAggregator
  next_jump::T
  end_time::T
  cur_rates::Vector{T}
  sum_rate::T
  rates::F1
  affects!::F2
  save_positions::Tuple{Bool,Bool}
end

@inline function (p::DirectVECJumpAggregation)(u,t,integrator) # condition
  p.next_jump==t
end

function (p::DirectVECJumpAggregation)(integrator) # affect!
  rng_val = rand() * p.sum_rate
  i       = searchsortedfirst(p.cur_rates, rng_val)
  @inbounds p.affects![i](integrator)
  p.sum_rate,ttnj = time_to_next_jump_ptr(integrator.u,integrator.p,integrator.t,p.rates,p.cur_rates)
  p.next_jump     = integrator.t + ttnj
  if p.next_jump < p.end_time
    add_tstop!(integrator,p.next_jump)
  end
  nothing
end

function (p::DirectVECJumpAggregation)(dj,u,t,integrator) # initialize
    sum_rate,next_jump = time_to_next_jump_ptr(u,integrator.p,t,p.rates,p.cur_rates)
    p.sum_rate         = sum_rate
    p.next_jump        = t + next_jump
    if p.next_jump < p.end_time
    add_tstop!(integrator,p.next_jump)
    end
    nothing
end

function time_to_next_jump_ptr(u,p,t,rates,cur_rates)
    @inbounds cur_rates[1] = rates[1](u,p,t)
    @inbounds sum_rate     = cur_rates[1]
    @inbounds for i in 2:length(cur_rates)
        new_rate     = rates[i](u,p,t)
        cur_rates[i] = new_rate + cur_rates[i-1]
        sum_rate    += new_rate
    end

    sum_rate,randexp()/sum_rate
end

@inline function aggregate(aggregator::DirectVEC,u,p,t,end_time,constant_jumps,save_positions)
    RateWrapper        = FunctionWrappers.FunctionWrapper{typeof(t),Tuple{typeof(u), typeof(p), typeof(t)}}
    rates              = [RateWrapper(c.rate) for c in constant_jumps]
    AffectWrapper      = FunctionWrappers.FunctionWrapper{Void,Tuple{Any}}
    affects!           = [AffectWrapper(x->(c.affect!(x);nothing)) for c in constant_jumps]
    # affects!           = [c.affect! for c in constant_jumps]
    cur_rates          = Vector{Float64}(length(rates))
    sum_rate,next_jump = time_to_next_jump_ptr(u,p,t,rates,cur_rates)
    DirectVECJumpAggregation(next_jump, end_time, cur_rates, sum_rate,
                                        rates, affects!, save_positions)
end

DiscreteCallback(c::DirectVECJumpAggregation) = DiscreteCallback(c,c,initialize=c,save_positions=c.save_positions)
