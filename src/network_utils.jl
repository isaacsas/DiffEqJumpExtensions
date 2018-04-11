using DataStructures

###############################################################################
# dep graph when MassActionJump uses pairs to represent (species,stoich)
###############################################################################

# map from species to Set of reactions depending on that species
function spec_to_dep_rxs_map(numspec, ma_jumps::MassActionJump)

    numrxs = length(ma_jumps.scaled_rates)

    # map from a species to reactions that depend on it
    spec_to_dep_rxs = [Set{Int}() for n = 1:numspec]
    for rx in 1:numrxs                    
        for (spec,stoch) in ma_jumps.reactant_stoch[rx]
            push!(spec_to_dep_rxs[spec], rx)
        end
    end

    spec_to_dep_rxs
end

# dependency graph is a map from a reaction to a vector of reactions
# that should depend on species it changes
function make_dependency_graph(numspec, ma_jumps::MassActionJump)

    numrxs          = length(ma_jumps.scaled_rates)
    spec_to_dep_rxs = spec_to_dep_rxs_map(numspec, ma_jumps)

    # create map from rx to reactions depending on it
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]
    for rx in 1:numrxs

        # rx changes spec, hence rxs depending on spec depend on rx
        for (spec,stoch) in ma_jumps.net_stoch[rx]
            for dependent_rx in spec_to_dep_rxs[spec]
                push!(dep_sets[rx], dependent_rx)
            end
        end
    end

    # convert to Vectors of Vectors
    dep_graph = Vector{Vector{Int}}(numrxs)
    for rx = 1:numrxs
        dep_graph[rx] = [dep for dep in dep_sets[rx]]
    end

    dep_graph
end


###############################################################################
# dep graph when MassActionJump uses 2xN matrix to represent (species,stoich)
###############################################################################


# # map from species to Set of reactions depending on that species
# function spec_to_dep_rxs_map(numspec, ma_jumps::MassActionJump)

#     numrxs = length(ma_jumps.scaled_rates)

#     # map from a species to reactions that depend on it
#     spec_to_dep_rxs = [Set{Int}() for n = 1:numspec]
#     for rx in 1:numrxs                    
#         for (spec,stoch) in ma_jumps.reactant_stoch[rx]
#             push!(spec_to_dep_rxs[spec], rx)
#         end
#     end

#     spec_to_dep_rxs
# end

# # dependency graph is a map from a reaction to a vector of reactions
# # that should depend on species it changes
# function make_dependency_graph(numspec, ma_jumps::MassActionJump)

#     numrxs          = length(ma_jumps.scaled_rates)
#     spec_to_dep_rxs = spec_to_dep_rxs_map(numspec, ma_jumps)

#     # create map from rx to reactions depending on it
#     dep_sets = [SortedSet{Int}() for n = 1:numrxs]
#     for rx in 1:numrxs

#         # rx changes spec, hence rxs depending on spec depend on rx
#         for (spec,stoch) in ma_jumps.net_stoch[rx]
#             for dependent_rx in spec_to_dep_rxs[spec]
#                 push!(dep_sets[rx], dependent_rx)
#             end
#         end
#     end

#     # convert to Vectors of Vectors
#     dep_graph = Vector{Vector{Int}}(numrxs)
#     for rx = 1:numrxs
#         dep_graph[rx] = [dep for dep in dep_sets[rx]]
#     end

#     dep_graph
# end
