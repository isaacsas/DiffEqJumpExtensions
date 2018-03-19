
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



@inbounds @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst,
                                        stochmat::AbstractVector{Pair{T,T}}) where T
    val = one(T)

    for specstoch in stochmat
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        for k = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

    return rateconst * val
end

@inline @inbounds @fastmath function executerx!(speciesvec::AbstractVector{T},
                                                net_stoch::AbstractVector{Pair{T,T}}) where T
    for specstoch in net_stoch
        speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end
