# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

# optimized version of DataStructure.jl priority queue
# swapped index from a dict to an array
# disabled checking for presence of keys when inserting!


import Base: <, <=, ==, length, isempty, start, next, done, delete!,
show, dump, empty!, getindex, setindex!, get, get!,
in, haskey, keys, merge, copy, cat,
push!, pop!, insert!,
union!, delete!, similar, sizehint!,
isequal, hash,
map, reverse,
first, last, eltype, getkey, values, sum,
merge, merge!, lt, Ordering, ForwardOrdering, Forward,
ReverseOrdering, Reverse, Lt,
isless,
union, intersect, symdiff, setdiff, issubset,
searchsortedfirst, searchsortedlast, in

import Base.Order: Forward, Ordering, lt

# Binary heap indexing
heapleft(i::Integer) = i << 1 #2i             # alt: i << 1
heapright(i::Integer) = (i << 1) + 1 #2i + 1        # alt: (i << 1) + 1
heapparent(i::Integer) = i >> 1 #div(i, 2)    # alt: i >> 1

function not_iterator_of_pairs(kv)
    return any(x->isempty(methodswith(typeof(kv), x, true)),
               [start, next, done]) ||
           any(x->!isa(x, Union{Tuple,Pair}), kv)
end

# ArrayPQ
# -------------

"""
    ArrayPQ(K, V, [ord])
Construct a new [`ArrayPQ`](@ref), with keys of type
`K` and values/priorites of type `V`.
If an order is not given, the priority queue is min-ordered using
the default comparison for `V`.
A `ArrayPQ` acts like a `Dict`, mapping integer values to their
priorities.
```jldoctest
julia> a = ArrayPQ([2, 4, 5],[2,3,1],Base.Order.Forward)
ArrayPQ{String,Int64,Base.Order.ForwardOrdering} with 3 entries:
  5 => 1
  4 => 3
  2 => 2
```
"""
mutable struct ArrayPQ{K,V,O<:Ordering} <: AbstractDict{K,V}
    # Binary heap of (element, priority) pairs.
    xs::Vector{Pair{K,V}}
    o::O

    # Map elements to their index in xs
    index::Vector{Int}

    function ArrayPQ{K,V,O}(o::O, itr) where {K,V,O<:Ordering}
        xs    = Vector{Pair{K,V}}(undef, length(itr))
        index = Vector{Int}(length(itr))
        @inbounds for (i,p) in enumerate(itr)
            xs[i]    = p
            index[i] = i
        end
        pq = new{K,V,O}(xs, o, index)

        # heapify
        for i in heapparent(length(pq.xs)):-1:1
            percolate_down!(pq, i)
        end

        pq
    end
end

function ArrayPQ(kv::Vector{T}, o::Ordering=Forward) where {K,V,T <: Pair{K,V}}
    try
        ArrayPQ{K,V,typeof(o)}(o, kv)
    catch e
        if not_iterator_of_pairs(kv)
            throw(ArgumentError("ArrayPQ(kv): kv needs to be an iterator of tuples or pairs"))
        else
            rethrow(e)
        end
    end
end
ArrayPQ(kv, o::Ordering=Forward) = ArrayPQ(o, kv)

length(pq::ArrayPQ) = length(pq.xs)
isempty(pq::ArrayPQ) = isempty(pq.xs)


"""
    peek(pq)
Return the lowest priority key from a priority queue without removing that
key from the queue.
"""
peek(pq::ArrayPQ) = pq.xs[1]

function percolate_down!(pq::ArrayPQ, i::Integer)
    @inbounds x = pq.xs[i]
    @inbounds while (l = heapleft(i)) <= length(pq)
        r = heapright(i)
        j = r > length(pq) || lt(pq.o, pq.xs[l].second, pq.xs[r].second) ? l : r
        if lt(pq.o, pq.xs[j].second, x.second)
            pq.index[pq.xs[j].first] = i
            pq.xs[i] = pq.xs[j]
            i = j
        else
            break
        end
    end
    @inbounds pq.index[x.first] = i
    @inbounds pq.xs[i] = x
end


function percolate_up!(pq::ArrayPQ, i::Integer)
    @inbounds x = pq.xs[i]
    @inbounds while i > 1
        j = heapparent(i)
        if lt(pq.o, x.second, pq.xs[j].second)
            pq.index[pq.xs[j].first] = i
            pq.xs[i] = pq.xs[j]
            i = j
        else
            break
        end
    end
    @inbounds pq.index[x.first] = i
    @inbounds pq.xs[i] = x
end

function getindex(pq::ArrayPQ{K,V}, key) where {K,V}
    @inbounds return pq.xs[pq.index[key]].second
end

# Change the priority of an existing element
# DOES NOT CHECK IF THE ELEMENT IS PRESENT!
function setindex!(pq::ArrayPQ{K, V}, value, key) where {K,V}
    @inbounds i = pq.index[key]
    @inbounds oldvalue = pq.xs[i].second
    @inbounds pq.xs[i] = Pair{K,V}(key, value)
    if lt(pq.o, oldvalue, value)
        percolate_down!(pq, i)
    else
        percolate_up!(pq, i)
    end
    value
end


# Unordered iteration through key value pairs in a ArrayPQ
start(pq::ArrayPQ) = start(pq.index)

done(pq::ArrayPQ, i) = done(pq.index, i)

function next(pq::ArrayPQ{K,V}, i) where {K,V}
    (k, idx), i = next(pq.index, i)
    return (pq.xs[idx], i)
end
