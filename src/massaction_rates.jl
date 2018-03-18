
@inbounds @fastmath function evalrxrate(speciesvec, rxidx, rateconstvec, stochmatvec)
    val      = rateconstvec[rxidx]  
    const stochmat = stochmatvec[rxidx]  
    for i = 1:2:length(stochmat)
        const specpop = speciesvec[stochmat[i]]
        val    *= specpop
        for k = 2:stochmat[i+1]
            val *= (specpop - k + 1)
        end
    end

    val
end

@inbounds @fastmath function executerx!(speciesvec, rxidx, net_stochvec)
    const net_stoch = net_stochvec[rxidx]
    for i = 1:2:length(net_stoch)
        speciesvec[net_stoch[i]] += net_stoch[i+1]
    end
    nothing 
end

