using FunctionWrappers

# when to switch aggregators from tuples for rates/affects to FunctionWrappers
const TUPLE_TO_FWRAPPER_CUTOFF = 10

function get_jump_info_tuples(constant_jumps)
    rates    = ((c.rate for c in constant_jumps)...)
    affects! = ((c.affect! for c in constant_jumps)...)
  
    return rates, affects!
  end
  
  function get_jump_info_fwrappers(u, p, t, constant_jumps)
    RateWrapper   = FunctionWrappers.FunctionWrapper{typeof(t),Tuple{typeof(u), typeof(p), typeof(t)}}
    rates         = [RateWrapper(c.rate) for c in constant_jumps]
    AffectWrapper = FunctionWrappers.FunctionWrapper{Void,Tuple{Any}}
    affects!      = [AffectWrapper(x->(c.affect!(x);nothing)) for c in constant_jumps]
  
    return rates, affects!
  end
