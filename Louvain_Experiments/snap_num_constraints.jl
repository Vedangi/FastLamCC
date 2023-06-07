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

function determine_class(graph)
    class = "--"
    classnum = 0
    if (graph in [ "loc-Brightkite";"loc-Gowalla"])
        class = "loc-social"
        classnum = 1
    elseif (graph in [ "com-LiveJournal";"com-Youtube";"soc-Epinions1";
        "soc-LiveJournal1";
        "soc-Slashdot0811";
        "soc-Slashdot0902"])
        class = "o-social"
        classnum = 2
    elseif (graph in [ "web-BerkStan";
        "web-Google";
        "web-NotreDame";
        "web-Stanford"; "wiki-topcats"])
        class = "web"
        classnum = 3
    elseif (graph in [     "email-Enron";
        "email-EuAll"; "wiki-Talk"])
        class = "comm"
        classnum = 4
    elseif (graph in [    "roadNet-CA";
        "roadNet-PA";
        "roadNet-TX"])
        class = "road"
        classnum = 5
    elseif (graph in [ "amazon0302";
        "amazon0312";
        "amazon0505";
        "amazon0601"; "com-Amazon" ])
        class = "prod"
        classnum = 6    
    elseif (graph in [ "com-DBLP";"ca-AstroPh";
        "ca-CondMat";
        "ca-GrQc";
        "ca-HepPh";
        "ca-HepTh"])
        class = "collab"
        classnum = 7 
    elseif (graph in [    "cit-HepPh";
        "cit-HepTh";
        "cit-Patents"])
        class = "cit"
        classnum = 8
    else
        println(" graph class not found")
    end
    return class,classnum
end

# Arrays to store the values
graph_results = Vector{Dict{String, Any}}()

 snaps = [
    "amazon0312";"amazon0302";
    "amazon0505";
    "amazon0601";
    "ca-AstroPh";
    "ca-CondMat";
    "ca-GrQc";
    "ca-HepPh";
    "ca-HepTh";
    "cit-HepPh";
    "cit-HepTh";
    "cit-Patents";
    "com-Amazon";
    "com-DBLP";
    "com-LiveJournal";
    "com-Youtube";
    "email-Enron";
    "email-EuAll";
    "loc-Brightkite";
    "loc-Gowalla";
    "roadNet-CA";
    "roadNet-PA";
    "roadNet-TX";
    "soc-Epinions1";
    "soc-LiveJournal1";
    "soc-Slashdot0811";
    "soc-Slashdot0902";
    "web-BerkStan";
    "web-Google";
    "web-NotreDame";
    "web-Stanford";
    "wiki-Talk"
    "wiki-topcats"
    ]

for i = 1:length(snaps)
    graph = snaps[i]
    F = matread("../datafiles/simple-snap/simple-$(graph).mat")
    A = Float64.(F["A"])
    n = size(A, 1)
    m = sum(A) / 2
    d = sum(diag(A))
    mx = maximum(A.nzval)``
    @assert(issymmetric(A))
    @assert(d == 0)

    if mx != 1
        Is, Js = findnz(A)
        A = sparse(Is, Js, 1, n, n)
        A = Float64.(A)
    end

    num_tri, num_ow = count_wedges_triangles(A)
    lp_constr = BigFloat(3 * binomial(BigInt(n), 3))

    class,classnum = determine_class(graph)

    graph_result = Dict("n" => n,"num_ow" => num_ow,"num_tri" => num_tri, "lp_constr" => lp_constr,"classnum" => classnum,"i" => i)
    file_path = "snap_constraints/$(graph)_dictionary_data.jld"

    # Serialize and save the dictionary to the file
    open(file_path, "w") do file
        serialize(file, graph_result)
    end
    
end


