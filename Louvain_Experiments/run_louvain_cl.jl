include("../include/lambda_louvain.jl")
include("../include/helpers.jl")
include("../include/faster_functions.jl")
include("../src/faster_functions_lamcc.jl")
include("../src/helpers_lamcc.jl")

"""
Run LambdaLouvain (lambda dependent) for one graph A.      

numtimes = number of times you run Louvain with different random node orderings.

maxits = number of iterations through nodes to allow
"""
function run_ll_cl(A,numtimes,maxits,lam)

    n = size(A,1)
    tic = time()
    c, lcc_obj = Many_Louvain(A,ones(n),lam,numtimes,maxits)
    timer = time()-tic

    Is, Js = findnz(triu(A))
    Elist = [Is Js]
    m = sum(A)/2

    mistakes = check_cl_obj_fastish(A,Elist,c,m,lam)
 
    return mistakes, timer
    
end
    