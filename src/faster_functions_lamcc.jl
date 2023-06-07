using Random
using SparseArrays
using LinearAlgebra
using MAT
include("../include/helpers.jl")
include("../include/faster_functions.jl")
include("../include/lp_functions.jl")
"""
Faster checking of LambdaCC correlation cluster objective.

Elist here is the actual edge list of a graph, don't double count.

It's not usually all pairs (unless this is a complete graph)
"""

function check_cl_obj_faster(Elist,c,m,Sbin,lam)
    n = length(c)
    pos_mis = 0
    for k = 1:size(Elist,1)
        if c[Elist[k,1]] != c[Elist[k,2]]
            pos_mis += 1
        end
    end

    neg_mis = pos_mis - m + Sbin
    return round(Int64,lam*neg_mis) + ((1-lam)*pos_mis)
end

"""
comparatively slow method to check LambdaCC objective 

"""
function cl_obj_slow(A,c,lam)
    n = size(A,1)
    mistakes = 0
    for i = 1:n
        for j = i+1:n
            if A[i,j] == 0 && c[i] == c[j]
                mistakes += lam
            end
            if A[i,j] == 1 && c[i] != c[j]
                mistakes += 1-lam
            end
        end
    end
    return mistakes
end


function many_pivot_cl(A,Elist,pivtimes ,lam)
    m = round(Int64,sum(A)/2)
    clus, Sbin = permutation_pivot_faster(A)
#     cc_obj = check_cc_obj_faster(Elist,clus,m,Sbin)
    cl_obj = check_cl_obj_faster(Elist,clus,m,Sbin,lam)
    # @assert(cc_obj == tst)

    # Pivot is fast, so we can run it multiple times
    for jj = 1:pivtimes
        clusnew, Sbin = permutation_pivot_faster(A)
#         objnew = check_cc_obj_faster(Elist,clusnew,m,Sbin)
        objnew = check_cl_obj_faster(Elist,clusnew,m,Sbin,lam)
        # @assert(objnew == tst)
        if objnew < cl_obj
            cl_obj = objnew
            clus = clusnew
        end
    end
    return cl_obj, clus
end


"""
Rounds a feasible LambdaSTC solution into a solution
    for a LambdaCC clustering problem.

Input:
    A = adjacency matrix

    Welist = node pairs to be flipped

Output:
    Approximation for CC.

    Step 1: Construct derived graph which is A with some flips
    Step 2: Apply pivot to the new graph
    Step 3: Compute CC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot
    method on the derived graph.

    Also returns the average time for a pivot step and the average
    quality solution
"""


function STCplus_to_CL_round_rand_faster(A,Eadd,Edel,pivtimes,lam)

    # Elist is list of edges for graph A
    n = size(A,1)
    tic = time()
    Anew = flip_graph(A,Eadd,Edel)
    setuptime = time()-tic
    # @show setuptime
    Is, Js = findnz(triu(A))
    ElistA = [Is Js]
    m = sum(A)/2

    NeighbsNew = ConstructAdj(Anew,n)

    # Apply pivot on this graph multiple times,
    # returning the best output

    clus = zeros(n)
    Sbin = permutation_pivot_fastest!(Anew,NeighbsNew,clus)
    cl_obj = check_cl_obj_faster(ElistA,clus,m,Sbin,lam)
    clusnew = zeros(n)
    objs = cl_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        # clusnew, Sbin = permutation_pivot_faster(Anew)
        # clusnew, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        Sbin = permutation_pivot_fastest!(Anew,NeighbsNew,clusnew)
        objnew = check_cl_obj_faster(ElistA,clusnew,m,Sbin,lam)
        
        # objnew2 = check_cc_obj(A,clusnew) # Check output on original graph
        # @show objnew, objnew2
        # @assert(objnew2 == objnew)

        objs += objnew
        if objnew < cl_obj
            cl_obj = objnew
            clus = clusnew
        end
    end
    totaltime = time()-tic
    avg_time = totaltime/pivtimes + setuptime
    avg_obj = objs/pivtimes
    return round(Int64,cl_obj), clus, avg_obj, avg_time
end