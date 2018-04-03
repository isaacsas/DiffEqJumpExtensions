
mutable struct DirectMAJumpAggregation{T,S,F1,F2,V} <: AbstractSSAJumpAggregator
    next_jump_time::T         # required field for SSAs, stores time of next jump
    next_jump::Int64          # required field for SSAs, stores index of next jump
    end_time::T               # required field for SSAs, stores time to stop simulation
    cur_rates::Vector{T}      # mass action rates followed by constant_jump rates
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool} # required field for SSAs, determines saving of output
    rng::V                    # required field for SSAs, stores random number generator
end

########### The following routines should be templates for all SSAs ###########

# condition for jump to occur
@inline function (p::DirectMAJumpAggregation)(u, t, integrator) 
    p.next_jump_time == t
end

# executing jump at the next jump time
function (p::DirectMAJumpAggregation)(integrator) 
    execute_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    generate_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    register_next_jump_time!(integrator, p, integrator.t)
    nothing
end

# setting up a new simulation
function (p::DirectMAJumpAggregation)(dj, u, t, integrator) # initialize
    initialize!(p, integrator, u, integrator.p, t)
    register_next_jump_time!(integrator, p, integrator.t)
    nothing
end


############################# Required Functions #############################

# creating the JumpAggregation structure
function aggregate(aggregator::DirectMassAction, u, p, t, end_time,
    constant_jumps, ma_jumps, save_positions; rng = Base.Random.GLOBAL_RNG)

    # handle constant jumps using function wrappers
    if constant_jumps != nothing
        rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    else
        rates    = []
        affects! = []
    end

    # current jump rates, allows mass action rates and constant jumps
    cur_rates = Vector{typeof(t)}(length(ma_jumps.scaled_rates) + length(rates))

    # setup the Jump Aggregation struct
    DirectMAJumpAggregation(zero(t), 0, end_time, cur_rates, zero(t), 
                            ma_jumps, rates, affects!, save_positions, rng)
end

# set up a new simulation and calculate the first jump / jump time
function initialize!(p::DirectMAJumpAggregation, integrator, u, params, t)
    generate_jumps!(p, integrator, u, integrator.p, t)    
end

# execute one jump, changing the system state
function execute_jumps!(p::DirectMAJumpAggregation, integrator, u, params, t)
    num_ma_rates = length(p.ma_jumps.scaled_rates)
    if p.next_jump <= num_ma_rates
        @inbounds executerx!(u, p.ma_jumps.net_stoch[p.next_jump])
    else
        idx = p.next_jump - num_ma_rates
        @inbounds p.affects![idx](integrator)
    end
end

# calculate the next jump / jump time
function generate_jumps!(p::DirectMAJumpAggregation, integrator, u, params, t)    
    # next jump time
    p.sum_rate, ttnj = time_to_next_jump_ma(p, u, params, t)
    @fastmath p.next_jump_time = t + ttnj    

    # next jump 
    rn = rand(p.rng) * p.sum_rate
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rn) 
end


######################## SSA specific helper routines ########################

@fastmath function time_to_next_jump_ma(p::DirectMAJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate  = zero(t)
    cur_rates = p.cur_rates

    # mass action rates
    majumps   = p.ma_jumps
    idx       = length(majumps.scaled_rates)
    @inbounds for i in 1:idx
        new_rate     = evalrxrate(u, majumps.scaled_rates[i], majumps.reactant_stoch[i])
        cur_rates[i] = new_rate + prev_rate
        prev_rate    = cur_rates[i]
    end

    # constant jump rates
    idx += 1
    @inbounds for i in 1:length(p.rates)
        new_rate        = p.rates[i](u, params, t)
        cur_rates[idx]  = new_rate + prev_rate
        prev_rate       = cur_rates[idx]
        idx            += 1
    end

    @inbounds sum_rate = cur_rates[end]
    sum_rate, randexp(p.rng) / sum_rate
end

