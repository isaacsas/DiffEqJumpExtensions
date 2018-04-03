

mutable struct DirectFunWrapJumpAgg{T,F1,F2} <: AbstractJumpAggregator
  next_jump_time::T
  next_jump::Int64
  end_time::T
  cur_rates::Vector{T}
  sum_rate::T
  rates::F1
  affects!::F2
  save_positions::Tuple{Bool,Bool}
end

@inline function (p::DirectFunWrapJumpAgg)(u,t,integrator) # condition
  p.next_jump_time == t
end

function (p::DirectFunWrapJumpAgg)(integrator) # affect!
  @inbounds p.affects![p.next_jump](integrator)
  update_samples(p, integrator, integrator.u, integrator.p, integrator.t)
  nothing
end

function (p::DirectFunWrapJumpAgg)(dj,u,t,integrator) # initialize
  update_samples(p, integrator, u, integrator.p, t)
  nothing
end

################ Direct specific implementation functions
function update_samples(p::DirectFunWrapJumpAgg, integrator, u, params, t)

  # next jump time
  p.sum_rate, ttnj = time_to_next_jump(u, params, t, p.rates, p.cur_rates)

  # add jump event to event queue
  p.next_jump_time = t + ttnj
  if p.next_jump_time < p.end_time
    add_tstop!(integrator, p.next_jump_time)
  end
 
  # next jump
  rn = rand() * p.sum_rate
  @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rn)
end


################# Tuple based time to next jump ##################
function time_to_next_jump(u, p, t, rates::Tuple, cur_rates) 
  @inbounds fill_cur_rates(u,p,t,cur_rates,1,rates...)
  @inbounds for i in 2:length(cur_rates) # cumsum
    cur_rates[i] = cur_rates[i] + cur_rates[i-1]
  end
  @inbounds sum_rate = cur_rates[end]
  sum_rate,randexp()/sum_rate
end

@inline function fill_cur_rates(u,p,t,cur_rates,idx,rate,rates...)
  @inbounds cur_rates[idx] = rate(u,p,t)
  idx += 1
  fill_cur_rates(u,p,t,cur_rates,idx,rates...)
end

@inline function fill_cur_rates(u,p,t,cur_rates,idx,rate)
  @inbounds cur_rates[idx] = rate(u,p,t)
  nothing
end

############ FunctionWrapper based time to next jump #############
function time_to_next_jump(u, p, t, rates::Vector{T}, cur_rates) where T
  prev_rate = zero(t)
  new_rate  = zero(t)
  @inbounds for i in 1:length(rates)
      new_rate     = rates[i](u, p, t)
      cur_rates[i] = new_rate + prev_rate
      prev_rate    = cur_rates[i]
  end
  @inbounds sum_rate = cur_rates[end]
  sum_rate,randexp()/sum_rate
end


function aggregate(aggregator::DirectFunWrappers, u, p, t, end_time, 
                           constant_jumps, ma_jumps, save_positions)
  
  # decide if representing rates/affects as tuples or function wrappers
  if length(constant_jumps) < TUPLE_TO_FWRAPPER_CUTOFF
    rates, affects! = get_jump_info_tuples(constant_jumps)  
  else
    rates, affects! = get_jump_info_fwrappers(u,p,t,constant_jumps)
  end

  cur_rates = zeros(typeof(t),length(rates))  
  DirectFunWrapJumpAgg(zero(t), 0, end_time, cur_rates, zero(t), rates, affects!, save_positions)
end

DiscreteCallback(c::DirectFunWrapJumpAgg) = DiscreteCallback(c, c, initialize=c, save_positions=c.save_positions)
