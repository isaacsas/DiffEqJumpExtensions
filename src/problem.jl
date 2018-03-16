DiffEqJump.JumpProblem(prob,aggregator::FRM,rn::DiffEqBiological.AbstractReactionNetwork;kwargs...) =
                 JumpProblem(prob,aggregator::FRM,rn.jumps...;kwargs...)

function DiffEqJump.JumpProblem(prob, aggregator::FRM, jumps::JumpSet;
    save_positions=typeof(prob) <: AbstractDiscreteProblem ? (false, true) : (true, true))
    
    ## Constant Rate Handling
    t, end_time, u = prob.tspan[1], prob.tspan[2], prob.u0
    if typeof(jumps.constant_jumps) <: Tuple{}
        disc = nothing
        constant_jump_callback = CallbackSet()
    else
        disc = aggregate(aggregator, u, prob.p, t, end_time, jumps.constant_jumps, save_positions)
        constant_jump_callback = DiscreteCallback(disc)
    end

    ## Variable Rate Handling
    if typeof(jumps.variable_jumps) <: Tuple{}
        new_prob = prob
        variable_jump_callback = CallbackSet()
    else
        new_prob = extend_problem(prob, jumps)
        variable_jump_callback = build_variable_callback(CallbackSet(), 0, jumps.variable_jumps...)
    end
    callbacks = CallbackSet(constant_jump_callback, variable_jump_callback)
    JumpProblem{typeof(new_prob),typeof(aggregator),typeof(callbacks),typeof(disc),typeof(jumps.variable_jumps)}(new_prob, aggregator, disc,
       callbacks,
       jumps.variable_jumps)
end


DiffEqJump.JumpProblem(prob,aggregator::FRMVEC,rn::DiffEqBiological.AbstractReactionNetwork;kwargs...) =
                 JumpProblem(prob,aggregator::FRMVEC,rn.jumps...;kwargs...)

function DiffEqJump.JumpProblem(prob, aggregator::FRMVEC, jumps::JumpSet;
    save_positions=typeof(prob) <: AbstractDiscreteProblem ? (false, true) : (true, true))
    
    ## Constant Rate Handling
    t, end_time, u = prob.tspan[1], prob.tspan[2], prob.u0
    if typeof(jumps.constant_jumps) <: Tuple{}
        disc = nothing
        constant_jump_callback = CallbackSet()
    else
        disc = aggregate(aggregator, u, prob.p, t, end_time, jumps.constant_jumps, save_positions)
        constant_jump_callback = DiscreteCallback(disc)
    end

    ## Variable Rate Handling
    if typeof(jumps.variable_jumps) <: Tuple{}
        new_prob = prob
        variable_jump_callback = CallbackSet()
    else
        new_prob = extend_problem(prob, jumps)
        variable_jump_callback = build_variable_callback(CallbackSet(), 0, jumps.variable_jumps...)
    end
    callbacks = CallbackSet(constant_jump_callback, variable_jump_callback)
    JumpProblem{typeof(new_prob),typeof(aggregator),typeof(callbacks),typeof(disc),typeof(jumps.variable_jumps)}(new_prob, aggregator, disc,
       callbacks,
       jumps.variable_jumps)
end

