using Random
using Gurobi
using SparseArrays
using LinearAlgebra
using MAT

# gurobi_env = Gurobi.Env()

"""
LambdaCC objective check

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

"""
Faster checking of LambdaCC objective
"""
function check_cl_obj_fastish(A::SparseArrays.SparseMatrixCSC,Elist,c,m,lam)
    n = length(c)
    pos_mis = 0
    for k = 1:size(Elist,1)
        if c[Elist[k,1]] != c[Elist[k,2]]
            pos_mis += 1
        end
    end

    clus_num = maximum(c)
    clus_sizes = zeros(Int64,clus_num)
    for i = 1:n
        clus_sizes[c[i]] += 1
    end

    neg_mis = (pos_mis - m)
    for k = 1:maximum(clus_num)
        S = clus_sizes[k]
        neg_mis += S*(S-1)/2
    end
    return round(Int64,(lam*neg_mis)) + ((1-lam)*pos_mis)
end





