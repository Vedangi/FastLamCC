# This is even faster code for many of the functions that rely on running pivot multiple times
# to either the original graph or some type of derived graph constructed based on
# a certain type of lower bound.
using Random
using SparseArrays
using LinearAlgebra
using MAT

function permutation_pivot_faster(A)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    Neighbs = ConstructAdj(A,n)
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end

"""
Fastest version of the pivot algorithm
"""
function permutation_pivot_fastest(A,Neighbs)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end

"""
Fastest version of the pivot algorithm
"""
function permutation_pivot_fastest!(A,Neighbs,c)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    # c = zeros(Int64,n)      # current cluster
    clusnum = 1

    for i = 1:n
        c[i] = 0
    end
    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return round(Int64,Sbin)
end
