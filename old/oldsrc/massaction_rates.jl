
###############################################################################
# Stochiometry for a given reaction is a vector of pairs mapping species id to
# stochiometric coefficient.
###############################################################################

@fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst,
                              stochmat::AbstractVector{Pair{S,V}})::typeof(rateconst) where {T,S,V}
    val = one(T)

    @inbounds for specstoch in stochmat
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        @inbounds for k = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

     rateconst * val
end

@inline @fastmath function executerx!(speciesvec::AbstractVector{T},
                                      net_stoch::AbstractVector{Pair{S,V}}) where {T,S,V}
    @inbounds for specstoch in net_stoch
        speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end


@inbounds function scalerates!(unscaled_rates, stochmat::Vector{Vector{Pair{S,T}}}) where {S,T}
    for i in eachindex(unscaled_rates)
        coef = one(T)
        for specstoch in stochmat[i]
            coef *= factorial(specstoch[2])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end


###############################################################################
# Stochiometry for a given reaction is a vector of 2xN matrices. 1st row is 
# species indices and second row is stoichiometry
###############################################################################

@fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst, stochmat::AbstractArray{S,2}) where {T,S}
    val = one(T)

    i = 1
    while i < length(stochmat)
        @inbounds specpop = speciesvec[stochmat[i]]
        val    *= specpop
        i      += 1
        @inbounds for k = 2:stochmat[i]
            specpop -= one(specpop)
            val     *= specpop
        end
        i += 1
    end

    return rateconst * val
end

@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, net_stoch::AbstractArray{S,2}) where {T,S}
    @inbounds for i = 1:2:length(net_stoch)
        speciesvec[net_stoch[i]] += net_stoch[i+1]
    end
    nothing
end

@inbounds function scalerates!(unscaled_rates, stochmat::Vector{Array{T,2}}) where {T}
    for i in eachindex(unscaled_rates)
        coef = one(T)
        stoch = stochmat[i]
        for k in 2:2:length(stoch)
            coef *= factorial(stoch[k])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end



# OLD!!!

# @inbounds @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst, stochmat::AbstractArray{T,2}) where T
#     val = one(T)

#     i = 1
#     while i < (length(stochmat))
#         specpop = speciesvec[stochmat[i]]
#         val    *= specpop
#         i      += 1
#         for k = 2:stochmat[i]
#             specpop -= one(specpop)
#             val     *= specpop
#         end
#         i += 1
#     end

#     return rateconst * val
# end

# @inline @inbounds @fastmath function executerx!(speciesvec::AbstractVector{T}, net_stoch::AbstractArray{T,2}) where T
#     for i = 1:2:length(net_stoch)
#         speciesvec[net_stoch[i]] += net_stoch[i+1]
#     end
#     nothing
# end



# @inbounds @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst,
#                                         stochmat::AbstractVector{Pair{T,T}}) where T
#     val = one(T)

#     for specstoch in stochmat
#         specpop = speciesvec[specstoch[1]]
#         val    *= specpop
#         for k = 2:specstoch[2]
#             specpop -= one(specpop)
#             val     *= specpop
#         end
#     end

#     return rateconst * val
# end

# @inline @inbounds @fastmath function executerx!(speciesvec::AbstractVector{T},
#                                                 net_stoch::AbstractVector{Pair{T,T}}) where T
#     for specstoch in net_stoch
#         speciesvec[specstoch[1]] += specstoch[2]
#     end
#     nothing
# end
