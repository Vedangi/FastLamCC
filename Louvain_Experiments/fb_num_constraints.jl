using MAT
using LinearAlgebra
using SparseArrays
using Statistics
using Combinatorics
using Serialization

include("../src/cfp_functions.jl")
include("../include/helpers.jl")
include("../include/faster_functions.jl")
using Plots

# Arrays to store the values
graph_results = Vector{Dict{String, Any}}()

fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")
 

for i = 1:100
    graph = fbs[i]
    F = matread("../datafiles/Facebook100/$(graph).mat")
    A = Float64.(F["A"])
    n = size(A, 1)
    m = sum(A) / 2
    d = sum(diag(A))
    mx = maximum(A.nzval)
    @assert(issymmetric(A))
    @assert(d == 0)

    if mx != 1
        Is, Js = findnz(A)
        A = sparse(Is, Js, 1, n, n)
        A = Float64.(A)
    end

    num_tri, num_ow = count_wedges_triangles(A)
    lp_constr = BigFloat(3 * binomial(BigInt(n), 3))


    graph_result = Dict("n" => n,"num_ow" => num_ow,"num_tri" => num_tri, "lp_constr" => lp_constr,"i" => i)
    file_path = "fb_constraints/$(graph)_dictionary_data.jld"

    # Serialize and save the dictionary to the file
    open(file_path, "w") do file
        serialize(file, graph_result)
    end
    
end


