
struct FRM <: AbstractAggregatorAlgorithm end
struct FRMVEC <: AbstractAggregatorAlgorithm end


# Mass action aggregators
abstract type AbstractMassActionAggregatorAlgorithm <: AbstractAggregatorAlgorithm end

mutable struct DirectMA{T,S} <: AbstractMassActionAggregatorAlgorithm 
    scaled_rates::T
    reactant_stoch::S
    net_stoch::S
end