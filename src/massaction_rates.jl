
@inbounds @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst, stochmat::AbstractArray{T,2}) where T
    val = one(T)

    i = 1
    while i < (length(stochmat))
        specpop = speciesvec[stochmat[i]]
        val    *= specpop
        i      += 1
        for k = 2:stochmat[i]
            specpop -= one(specpop)
            val     *= specpop
        end
        i += 1
    end

    return rateconst * val
end

@inline @inbounds @fastmath function executerx!(speciesvec::AbstractVector{T}, net_stoch::AbstractArray{T,2}) where T
    for i = 1:2:length(net_stoch)
        speciesvec[net_stoch[i]] += net_stoch[i+1]
    end
    nothing
end
