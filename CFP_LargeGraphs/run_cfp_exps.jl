using MAT
using LinearAlgebra
using SparseArrays
using Statistics

include("../src/cfp_functions.jl")
include("../src/helpers_lamcc.jl")
include("../src/faster_functions_lamcc.jl")
include("../include/helpers.jl")
include("../include/faster_functions.jl")
include("../include/lp_functions.jl")

snaps = ["amazon0302";
    "amazon0312";
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
    "wiki-Talk";
    "wiki-topcats"
    ]


# fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")
# graph = fbs[8]
# F = matread("../datafiles/Facebook100/$graph.mat")

graph = snaps[29]
F = matread("../datafiles/simple-snap/simple-$(graph).mat")


A = Float64.(F["A"])
n = size(A,1)
m = sum(A)/2
d = sum(diag(A))
mx = maximum(A.nzval)
@assert(issymmetric(A))
@assert(d == 0)
if mx != 1
    Is, Js = findnz(A)
    A = sparse(Is,Js,1,n,n)
    A = Float64.(A)
end

println("$graph \t $n \t $m")

num_runs = 15
lam_values = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99]
for lam in lam_values
    # Get lower bound , edges to delete and add
    tic = time()
    Edel, Eadd, lb = coverLabel(A,lam)
    lb_time = time()-tic #lower bound time
    
    cl_ratio_results = zeros(15)
    cl_ub_results = zeros(15)
    cl_ubt_results = zeros(15)
    cl_tott_results = zeros(15)
    
    graph_stats = [n;m;lam]
 
    pivtimes = 100
    
    for j in 1:num_runs

        # Randomized rounding
        tic = time()
        ub_rand_best, c_rand, ub_rand_avg, time_rand_avg = STCplus_to_CL_round_rand_faster(A,Eadd, Edel,pivtimes,lam)
        time_rand = time()-tic
        
         
        cl_ub_results[j] = ub_rand_best # cfp upper bound
        cl_ubt_results[j] = time_rand #upper bound time
        
        cl_ratio_results[j] =ub_rand_best/lb   # cfp ratio (upper bound/lower bound)
        cl_tott_results[j] = lb_time+time_rand # cfp runtime(lower bound time + upper bound time)
     
    end

    #Store Facebook graph results in "fb_cfp/$(graph)_fb_cfp_lam_$(lam).mat"
    #Store SNAP graph results in "snap_cfp/$(graph)_fb_cfp_lam_$(lam).mat"
    matwrite("snap_cfp/$(graph)_fb_cfp_lam_$(lam).mat", Dict("graph_stats"=>graph_stats,"lb"=>lb,"lb_time"=>lb_time,"cl_ratio_results"=>cl_ratio_results,"cl_ub_results"=>cl_ub_results,"cl_ubt_results"=>cl_ubt_results,"cl_tott_results"=>cl_tott_results))

end