using Random
using SparseArrays
using LinearAlgebra
using MAT

"""
Given a set of edges to add and delete, create a
new graph with all these edges flipped
"""
function flip_graph(A,Eadd,Edel)
    n = size(A,1)
    Anew = copy(A)
    for k = 1:size(Edel,1)
        i = Edel[k][1]
        j = Edel[k][2]
        Anew[i,j] = 0
        Anew[j,i] = 0
    end
    dropzeros!(Anew)
    II,JJ = findnz(Anew)
    for k = 1:size(Eadd,1)
        push!(II,Eadd[k][1])
        push!(JJ,Eadd[k][2])
        push!(II,Eadd[k][2])
        push!(JJ,Eadd[k][1])
    end
    Anew = sparse(II,JJ,1,n,n)

    return Anew
end

function STCplus_to_CL_round_rand(A,Eadd, Edel,pivtimes,lam)

    tic = time()
    Anew = flip_graph(A,Eadd, Edel)
    setuptime = time()-tic

    # Apply pivot on this graph multiple times,
    # returning the best output
    clus = permutation_pivot(Anew)
#     cl_obj = check_cl_obj(A,clus,lam)
    cl_obj = cl_obj_slow(A,clus,lam)
    objs = cl_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        clusnew = permutation_pivot(Anew)
#         objnew = check_cl_obj(A,clus,lam) # Check output on original graph
        objnew = cl_obj_slow(A,clus,lam)
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



"""
Standard rounding technique 

Rounds a lower bound for LambdaSTC into a feasible solution
    for a LambdaCC clustering.

This works for:
    * rounding the LambdaCC LP relaxation (3-approx)
    * rounding the LambdaSTC LP relaxation (6-approx)
    * rounding a feasible solution to the LambdaSTC ILP (6-approx)

Input:
    A = adjacency matrix

    Elist = linear ordering of node-pairs

    soln = feasible solution to the CC LP relaxation
            or the LambdaSTC LP relaxation
            or the LambdaSTC ILP

Output:
    Approximation for LambdaCC.

    Step 1: Construct derived graph where edges are pairs with x[e] < 1/3
    Step 2: Apply pivot to the new graph
    Step 3: Compute LambdaCC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot 
    method on the derived graph.
""" 

function CL_round_LP(A,Elist,soln,pivtimes,lam)
    n = size(A,1)

    # Elist is a linearization of all node pairs in A,
    # not just the edges in A.

    # Find all the node pairs Elist[k] that have soln[k] < 1/2
    tic = time()
    tol = 1e-8
    thr = 1/3
    keep = findall(x->x<thr-tol,soln)

    # These all become edges in a new derived graph Anew
    Inew = Elist[keep,1]
    Jnew = Elist[keep,2]
    Anew = sparse(Inew,Jnew,ones(length(Jnew)),n,n)
    Anew = Anew + Anew'

    Is, Js = findnz(triu(A))
    ElistA = [Is Js]
    m = sum(A)/2

    NeighbsNew = ConstructAdj(Anew,n)
    setuptime = time()-tic
    # Apply pivot on this graph multiple times,
    # returning the best output

    clus, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
    cl_obj = check_cl_obj_faster(ElistA,clus,m,Sbin,lam)

    objs = cl_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times

        clusnew,Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        objnew = check_cl_obj_faster(ElistA,clusnew,m,Sbin,lam)

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