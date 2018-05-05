# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

# simple array based priority queue, interface similar to 
# DataStructures.jl prioirty queue, but simplified
# disabled checking for presence of keys when inserting

# Binary heap indexing
heapleft(i::Integer)   = i << 1        #2i            
heapright(i::Integer)  = (i << 1) + 1 #2i + 1 
heapparent(i::Integer) = i >> 1      #div(i, 2)   


mutable struct SimplePQ{K,V} 
    # Binary heap of (element, priority) pairs.
    xs::Vector{Pair{K,V}}

    # Map elements to their index in xs
    index::Vector{Int}

    function SimplePQ{K,V}(itr) where {K,V}
        xs    = Vector{Pair{K,V}}(undef, length(itr))
        index = collect(1:length(itr))
        @inbounds for (i,p) in enumerate(itr)
            xs[i] = p
        end
        pq = new{K,V}(xs, index)

        # heapify
        for i in heapparent(length(pq.xs)):-1:1
            heapify(pq, i)
        end

        pq
    end
end

SimplePQ(kv::Vector{T}) where {K,V,T <: Pair{K,V}} = SimplePQ{K,V}(kv)   
length(pq::SimplePQ)  = length(pq.xs)
isempty(pq::SimplePQ) = isempty(pq.xs)
peek(pq::SimplePQ)    = (@inbounds return pq.xs[1])
@inline function swap(pq::SimplePQ{K,V}, i::Int, j::Int) where {K,V}
    @inbounds pq.xs[i], pq.xs[j] = pq.xs[j], pq.xs[i]    
    @inbounds ki = pq.xs[j].first
    @inbounds kj = pq.xs[i].first
    @inbounds pq.index[ki], pq.index[kj] = pq.index[kj], pq.index[ki]
end

@inbounds function update_pq_at(pq::SimplePQ{K,V}, i::Int) where {K,V}
    pidx = heapparent(i)
    xs   = pq.xs
    if (pidx > 0) && xs[i].second < xs[pidx].second
        swap(pq, i, pidx)
        update_pq_at(pq, pidx)
    else
        ridx = heapright(i)
        lidx = heapleft(i)
        len  = length(pq)
        if ridx <= len
            if (lidx <= len) && (xs[lidx].second < xs[ridx].second)
                ridx = lidx
            end
        else
            if lidx > len
                return
            end
            ridx = lidx
        end

        if xs[ridx].second < xs[i].second
            swap(pq, i, ridx)
            update_pq_at(pq, ridx)
        end
    end
end

@inbounds function heapify(pq::SimplePQ{K,V}, i::Int) where {K,V}
    l   = heapleft(i)
    r   = heapright(i)
    len = length(pq)
    xs  = pq.xs
    
    smallest = ((l <= len) && (xs[l].second < xs[i].second)) ? l : i

    if (r <= len) && (xs[r].second < xs[smallest].second)
        smallest = r
    end

    if smallest != i
        swap(pq, i, smallest)
        heapify(pq, smallest)        
    end
end


function getindex(pq::SimplePQ{K,V}, key) where {K,V}
    @inbounds return pq.xs[pq.index[key]].second
end

# Change the priority of an existing element
function setindex!(pq::SimplePQ{K, V}, value, key) where {K,V}
    @inbounds i        = pq.index[key]
    @inbounds pq.xs[i] = Pair{K,V}(key, value)
    update_pq_at(pq, i)
    value
end


# Unordered iteration through key value pairs in a SimplePQ
start(pq::SimplePQ)   = start(pq.index)
done(pq::SimplePQ, i) = done(pq.index, i)
function next(pq::SimplePQ{K,V}, i) where {K,V}
    (k, idx), i = next(pq.index, i)
    return (pq.xs[idx], i)
end
