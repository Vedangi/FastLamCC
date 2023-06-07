# include("../src/coverflip_functions.jl")
include("../src/cfp_functions.jl")
include("../src/helpers_lamcc.jl")
include("../src/faster_functions_lamcc.jl")
include("../src/lp_relaxations_lamcc.jl")
include("../src/rounding.jl")
include("../include/helpers.jl")
include("../include/faster_functions.jl")
include("../include/lp_functions.jl")


using Gurobi
using SparseArrays
using LinearAlgebra
using MAT

GC.gc()

gurobi_env = Gurobi.Env()

tinytestgraphs = [
    "KarateA";
    "dolphinsA";
    "lesmisA";
    "polbooksA";
    "adjnounA";
    "footballA";
]

smallgraphs = [
    "Harvard500A";
    "Erdos991A"; 
    "celegansneuralA";
    "Netscience";
    "celegansmetabolicA";
    "RogetA";
    "SmaGriA";
    "emailA";
    "polblogsA";
    "standard_ca-GrQc";
    "standard_caHepThA";
    "standard_EmailEnronA";
    "standard_condmat2005A";
    "standard_ca-AstroPhA";
    "standard_loc-Brightkite";
    # "Caltech36";
    # "Reed98";
    # "Simmons81";
    # "Haverford76";
    # "Swarthmore42";
    # "Amherst41";
    # "Bowdoin47";
    # "Rice31";
    # "Lehigh96";
]
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

# Pick a graph from the list


# fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")
# graph = fbs[56]
# graph = snaps[3]
graph = smallgraphs[6]

tl = 7200           # time limit
pivtimes = 100      # number of times to run pivot in rounding step
   
F = matread("../data/smallgraphs/$graph.mat")
# F = matread("../datafiles/simple-snap/simple-$graph.mat")
# F = matread("../datafiles/Facebook100/$graph.mat")

A = F["A"]
n = size(A,1)
# A = sparse(A)
n = size(A,1)
m = sum(A)/2
d = sum(diag(A))
mx = maximum(A.nzval)
@assert(issymmetric(A))
@assert(d == 0)
if mx != 1
    Is, Js = findnz(A)
    A = sparse(Is,Js,1,n,n)
end
num_runs = 15
lam_values = [0.4,0.55,0.75,0.85,0.9,0.95,0.99] #

#Calculating CFP 
for lam in lam_values  
    
    # Get lower bound and edges to be deleted and added
    tic = time()
    Edel, Eadd, lb = coverLabel(A,lam)
    lb_time = time()-tic #cfp lower bound time
    
    cfp_ratio_results = zeros(15)
    cfp_ub_results = zeros(15)
    cfp_ubt_results = zeros(15)
    cfp_tott_results = zeros(15)

    
    graph_stats = [n;m;lam]
 
    pivtimes = 100
    
    for j in 1:num_runs

        # Randomized rounding
        tic = time()
        ub_rand_best, c_rand, ub_rand_avg, time_rand_avg = STCplus_to_CL_round_rand_faster(A,Eadd, Edel,pivtimes,lam)
        time_rand = time()-tic 
         
        cfp_ub_results[j] = ub_rand_best #cfp upper bound
        cfp_ubt_results[j] = time_rand #cfp upper bound time
        
        cfp_ratio_results[j] =ub_rand_best/lb        #cfp ratio (upper bound/lower bound)
        cfp_tott_results[j] = lb_time+time_rand #cfp total time(lower bound time + upper bound time)
     
 
    end
    #Save CFP results
    matwrite("cfp_smallgraph_results/$(graph)_cfp_combined_lam_$(lam).mat", Dict("graph_stats"=>graph_stats,"lambda"=> lam, "lb"=>lb,"lb_time"=>lb_time,"cfp_ratio_results"=>cfp_ratio_results,"cfp_ub_results"=>cfp_ub_results,"cfp_ubt_results"=>cfp_ubt_results,"cfp_tott_results"=>cfp_tott_results))

end

#For LambdaCC LP and LambdaSTC LP relaxations
println("LP calculations")
for lam in lam_values
    
    time_limit = 7200
    pivtimes = 100 
    
    println("$graph \t lambda- $lam \t $n \t $m")
    graph_stats = [n;m;lam]
    
    if m < 500000
        checkcc = true
    else
        checkcc = false
    end
    
  
    lp_lb_results = zeros(15)
    lp_lbt_results = zeros(15)   
    clp_lb_results = zeros(15)
    clp_lbt_results = zeros(15)
    
#     lp_std_ratio_results = zeros(15)
    lp_std_ub_results = zeros(15)
    lp_std_tott_results = zeros(15)
    lp_std_ubt_results = zeros(15)
    
#     lp_new_ratio_results= zeros(15)
    lp_new_ub_results = zeros(15)
    lp_new_tott_results = zeros(15)
    lp_new_ubt_results = zeros(15)
              
#     clp_ratio_results = zeros(15)
    clp_new_ub_results = zeros(15)      
    clp_new_tott_results = zeros(15)
    clp_new_ubt_results = zeros(15)

    # clp_std_ratio_results = zeros(15)
    clp_std_ub_results = zeros(15)      
    clp_std_tott_results = zeros(15)
    clp_std_ubt_results = zeros(15)
    
    outflag = false
    LP = true
    verb = true
    tic = time()
    status_stclp, lb_stclp, stclp_soln, lb_time_stclp, status_cclp, lb_cclp, Elist, cclp_soln = LazyCLandcheap_combined(A;outputflag = outflag, verbose = verb,LP = LP,timelimit = tl,lam=lam)
    lb_time_cclp = time()-tic
    println("Status cclp is for $(status_cclp) for lambda $(lam)")
        # Does it satisfy CL constraints?
    if checkcc
        iscclp = is_CL_feasible(A,Elist,stclp_soln)
    else
        iscclp = -1
    end
 
    
    
    for j in 1:num_runs
        
          # round LambdaCC LP standard technique
        tic = time()
        cclp_round_std,c_std = CL_round_LP(A,Elist,cclp_soln,pivtimes,lam)
        time_lp_rd_std = time()-tic
        
          # round LambdaCC LP using new twchnique
        tic = time()
        cclp_round_new,c_new = CL_round_LP_new(A,Elist,cclp_soln,pivtimes,lam)
        time_lp_rd_new = time()-tic
        
        # round LambdaSTC LP  standard procedure
        tic = time()
        ub_stclp_std, c_std, avg_stclp_std, avgtime_stclp_std = CL_round_LP(A,Elist,stclp_soln,pivtimes,lam)
        time_round_stclp_std = time()-tic

        # round LambdaSTC LP new technique
        tic = time()
        ub_stclp_new, c_new, avg_stclp_new, avgtime_stclp_new = CL_round_LP_new(A,Elist,stclp_soln,pivtimes,lam)
        time_round_stclp_new = time()-tic
        
        #-------------------------------------------------------------------------------------------------------
        lp_lb_results[j] = lb_cclp  #LambdaCC LP lower bound
        lp_lbt_results[j] = lb_time_cclp #LambdaCC LP lower bound time
        
        clp_lb_results[j] = lb_stclp  #LambdaSTC LP lower bound
        clp_lbt_results[j] = lb_time_stclp #LambdaSTC LP lower bound time
        
#         lp_std_ratio_results[j] = cclp_round_std/lb_cclp
        lp_std_ub_results[j] = cclp_round_std #LambdaCC LP standard rounding upper bound
        lp_std_tott_results[j] = lb_time_cclp+time_lp_rd_std #LambdaCC LP standard rounding total time
        lp_std_ubt_results[j] = time_lp_rd_std #LambdaCC LP standard rounding upper bound time
        
#         lp_new_ratio_results[j] = cclp_round_new/lb_cclp
        lp_new_ub_results[j] = cclp_round_new #LambdaCC LP new rounding upper bound
        lp_new_tott_results[j] = lb_time_cclp+time_lp_rd_new #LambdaCC LP new rounding total time
        lp_new_ubt_results[j] = time_lp_rd_new #LambdaCC LP new rounding upper bound time
              
#        LambdaSTC LP standard rounding results
        clp_std_ub_results[j] = ub_stclp_std   #LambdaSTC LP standard rounding upper bound
        clp_std_tott_results[j] = lb_time_stclp+time_round_stclp_std #LambdaSTC LP standard total time
        clp_std_ubt_results[j] = time_round_stclp_std #LambdaSTC LP standard rounding upper bound time

        # LambdaSTC LP new rounding method results
        clp_new_ub_results[j] = ub_stclp_new  #LambdaSTC LP new rounding upper bound
        clp_new_tott_results[j] = lb_time_stclp+time_round_stclp_new #LambdaSTC LP new rounding total time
        clp_new_ubt_results[j] = time_round_stclp_new #LambdaSTC LP new rounding upper bound time
    
    end
   
    
    matwrite("lp_smallgraph_results/$(graph)_lp_combined_lam_$(lam).mat", Dict("pivtimes"=>pivtimes,"iscclp"=>iscclp,"status_stclp"=>status_stclp,"status_cclp"=>status_cclp,"lp_lb_results"=>lp_lb_results,"lp_lbt_results"=>lp_lbt_results,"clp_lb_results"=>clp_lb_results,"clp_lbt_results"=>clp_lbt_results,"lp_std_ub_results"=>lp_std_ub_results,"lp_std_tott_results"=>lp_std_tott_results,"lp_std_ubt_results"=>lp_std_ubt_results,"lp_new_ub_results"=>lp_new_ub_results,"lp_new_tott_results"=>lp_new_tott_results,"lp_new_ubt_results"=>lp_new_ubt_results,"clp_std_ub_results"=>clp_std_ub_results,"clp_std_tott_results"=>clp_std_tott_results,"clp_std_ubt_results"=>clp_std_ubt_results,"clp_new_ub_results"=>clp_new_ub_results,"clp_new_tott_results"=>clp_new_tott_results,"clp_new_ubt_results"=>clp_new_ubt_results))
    # matwrite("lp_snap_results/$(graph)_lp_combined_lam_$(lam).mat", Dict("pivtimes"=>pivtimes,"iscclp"=>iscclp,"status_stclp"=>status_stclp,"status_cclp"=>status_cclp,"lp_lb_results"=>lp_lb_results,"lp_lbt_results"=>lp_lbt_results,"clp_lb_results"=>clp_lb_results,"clp_lbt_results"=>clp_lbt_results,"lp_std_ub_results"=>lp_std_ub_results,"lp_std_tott_results"=>lp_std_tott_results,"lp_std_ubt_results"=>lp_std_ubt_results,"lp_new_ub_results"=>lp_new_ub_results,"lp_new_tott_results"=>lp_new_tott_results,"lp_new_ubt_results"=>lp_new_ubt_results,"clp_std_ub_results"=>clp_std_ub_results,"clp_std_tott_results"=>clp_std_tott_results,"clp_std_ubt_results"=>clp_std_ubt_results,"clp_new_ub_results"=>clp_new_ub_results,"clp_new_tott_results"=>clp_new_tott_results,"clp_new_ubt_results"=>clp_new_ubt_results))
    # matwrite("lp_fb_results/$(graph)_lp_combined_lam_$(lam).mat", Dict("pivtimes"=>pivtimes,"iscclp"=>iscclp,"status_stclp"=>status_stclp,"status_cclp"=>status_cclp,"lp_lb_results"=>lp_lb_results,"lp_lbt_results"=>lp_lbt_results,"clp_lb_results"=>clp_lb_results,"clp_lbt_results"=>clp_lbt_results,"lp_std_ub_results"=>lp_std_ub_results,"lp_std_tott_results"=>lp_std_tott_results,"lp_std_ubt_results"=>lp_std_ubt_results,"lp_new_ub_results"=>lp_new_ub_results,"lp_new_tott_results"=>lp_new_tott_results,"lp_new_ubt_results"=>lp_new_ubt_results,"clp_std_ub_results"=>clp_std_ub_results,"clp_std_tott_results"=>clp_std_tott_results,"clp_std_ubt_results"=>clp_std_ubt_results,"clp_new_ub_results"=>clp_new_ub_results,"clp_new_tott_results"=>clp_new_tott_results,"clp_new_ubt_results"=>clp_new_ubt_results))


end