using Gurobi
using SparseArrays
using LinearAlgebra
"""
New type of rounding technique 

Rounds a lower bound for LambdaSTC into a feasible solution
    for a LambdaCC clustering.

This works for:
    * rounding the LambdaCC LP relaxation
    * rounding the LambdaSTC LP relaxation for lambda>=0.5 (7 - 2/lambda approx)
    * rounding the LambdaSTC LP relaxation for lambda < 0.5 ((1+lambda)/lambda - approx)

Input:
    A = adjacency matrix

    Elist = linear ordering of node-pairs

    soln = feasible solution to the CC LP relaxation
            or the LambdaSTC LP relaxation
            or the LambdaSTC ILP

Output:
    Approximation for LambdaCC.

    Step 1: Construct derived graph based on the new rounding technique
    Step 2: Apply pivot to the new graph
    Step 3: Compute LambdaCC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot 
    method on the derived graph.
"""

function CL_round_LP_new(A,Elist,soln,pivtimes,lam)
    n = size(A,1)
   
    # Elist is a linearization of all node pairs in A,
    # not just the edges in A.

    tic = time()
    tol = 1e-8
    thr = (2*lam)/((7*lam)-2)
    
    # Select the elements of soln that correspond to edges/non-edges in A according to the rounding technique
    keep = Int64[]

    if (lam < 0.5)
        thr = (lam)/(lam+1)
        for k in 1:size(Elist, 1)
            i = Elist[k,1]
            j = Elist[k,2]
            if i < j
                if ((A[i, j] != 0) || (soln[k] < (thr - tol)))  # Find the indices of the non-edge pairs that are less than thr-tol
                    push!(keep, k)
                end
            end
        end

    else
        for k in 1:size(Elist, 1)
            i = Elist[k,1]
            j = Elist[k,2]
            if i < j && soln[k] < (thr - tol)  # Find the indices of the edge pairs that are less than thr-tol
                if A[i, j] != 0
                    push!(keep, k)
                end
            end
        end
    end

   
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

