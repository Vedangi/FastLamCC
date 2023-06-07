using MAT
using Statistics
using SparseArrays

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
    # "Simmons81"; facebook
    # "Haverford76";
    # "Swarthmore42";
    # "Amherst41";
    # "Bowdoin47"; facebook grph
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
fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")

# graph = fbs[56]
# graph = snaps[3]

graph = smallgraphs[6]

# F = matread("../datafiles/Facebook100/$graph.mat")
F = matread("../data/smallgraphs/$graph.mat")
# F = matread("../datafiles/simple-snap/simple-$graph.mat")

A = F["A"]
n = size(A,1)
m = round(Int64,sum(A)/2)
dig = 2
dig2 = 2

println("n($n) and m($m) for $graph ")
## collect
function prun(runtime,dig)
    if runtime > 1
        rt = round(runtime,digits = 1)
    else
        rt = round(runtime, sigdigits = dig)
    end
    return rt
end

for lam in [0.9,0.95,0.99] #0.4,0.55,0.75,0.85,

   
    # Print Combinatorial cfp results

    M = matread("cfp_smallgraph_results/$(graph)_cfp_combined_lam_$(lam).mat")
    #  M = matread("cfp_snap_results/$(graph)_cfp_combined_lam_$(lam).mat")
    # M = matread("cfp_fb_results/$(graph)_cfp_combined_lam_$(lam).mat")
   
    
    cfp_ub_mean = mean(M["cfp_ub_results"]) #cfp upper bound mean
    cfp_ub_std = std(M["cfp_ub_results"]) #cfp upper bound std
    cfp_lb = M["lb"] #cfp lower bound 
    cfp_lbt = M["lb_time"] #cfp lower bound  time
    cfp_ratio_mean = mean(M["cfp_ratio_results"]) #cfp ratio mean
    cfp_ratio_std = std(M["cfp_ratio_results"]) #cfp ratio std
      
    cfp_ubt_mean = mean(M["cfp_ubt_results"]) #cfp upper bound time mean
    cfp_ubt_std = std(M["cfp_ubt_results"]) #cfp upper bound time std

    cfp_tott_mean = mean(M["cfp_tott_results"]) #cfp total time mean
    cfp_tott_std = std(M["cfp_tott_results"])#cfp total time std


#     Printing rounded values
    
    lb_cfp = round(Int64,cfp_lb)
    ub_cfp_mean_round = round(Int64,cfp_ub_mean)
    ub_cfp_std_round = round(Int64,cfp_ub_std)

    rat_cfp_mean_round = round(cfp_ratio_mean,sigdigits = dig2)
    rat_cfp_std_round = round(cfp_ratio_std,sigdigits = dig2)
    
    run_cfp_mean_round = prun((cfp_tott_mean),dig)
    run_cfp_std_round = prun((cfp_tott_std),dig)

# LambdaSTC LP 
    oom = "--"
    lb_cheaplp = oom
    ub_cheaplp_rand_mean = oom
    ub_cheaplp_rand_std = oom
    rat_cheaplp_rand = oom
    run_cheaplp_rand = oom
    status_cheaplp = 1.0
    iscc_cheaplp = 0.0
    status_lp = 1.0
    run_lp = 0.0
    lb_lp = 0.0
    ub_lp = 0.0
    rat_lp = 0.0
    try
        M = matread("lp_smallgraph_results/$(graph)_lp_combined_lam_$(lam).mat")
        # M = matread("lp_snap_results/$(graph)_lp_combined_lam_$(lam).mat")
        # M = matread("lp_fb_results/$(graph)_lp_combined_lam_$(lam).mat")
        

#        
        lb_cheaplp = round(mean(M["clp_lb_results"]),digits = 1) #LambdaSTC LP lower bound
        
        #standard rounding procedure using x_{ij} < 1/3
        ub_cheaplp_rand_std_mean = round(Int64, mean(M["clp_std_ub_results"])) #standard rounding LambdaSTC LP upper bound mean
        ub_cheaplp_rand_std_std = round(Int64,std(M["clp_std_ub_results"])) #standard rounding LambdaSTC LP upper bound std
        rat_cheaplp_rand_std_mean = round(mean(M["clp_std_ub_results"] ./ M["clp_lb_results"] ),sigdigits=dig2) #standard rounding LambdaSTC LP ratio mean
        rat_cheaplp_rand_std_std = round(std(M["clp_std_ub_results"] ./ M["clp_lb_results"] ),sigdigits = dig2) #standard rounding LambdaSTC LP ration std

        run_cheaplp_rand_std_mean = prun(mean(M["clp_std_tott_results"]),dig) #standard rounding LambdaSTC LP total runtime mean
        run_cheaplp_rand_std_std = prun(std(M["clp_std_tott_results"]),dig) #standard rounding LambdaSTC LP total runtime mean
        
        #new rounding cheaplp using x_{ij} < 2lambda/(7lambda - 2 )
        ub_cheaplp_rand_new_mean = round(Int64, mean(M["clp_new_ub_results"])) #new rounding LambdaSTC LP upper bound mean
        ub_cheaplp_rand_new_std = round(Int64,std(M["clp_new_ub_results"])) #new rounding LambdaSTC LP upper bound std
        rat_cheaplp_rand_new_mean = round(mean(M["clp_new_ub_results"] ./ M["clp_lb_results"] ),sigdigits=dig2) #new rounding LambdaSTC LP ratio mean
        rat_cheaplp_rand_new_std = round(std(M["clp_new_ub_results"] ./ M["clp_lb_results"] ),sigdigits = dig2) #new rounding LambdaSTC LP ratio std

        run_cheaplp_rand_new_mean = prun(mean(M["clp_new_tott_results"]),dig) #new rounding LambdaSTC LP total runtime mean
        run_cheaplp_rand_new_std = prun(std(M["clp_new_tott_results"]),dig) #new rounding LambdaSTC LP total runtime std


        
        status_cheaplp = M["status_stclp"] # status--if false, we didn't converge
        iscc_cheaplp = M["iscclp"]   # if true, this is the LambdaCC solution too!
        
#       LambdaCC LP
        
        lb_lp = round(mean(M["lp_lb_results"]),digits = 1) #LambdaCC LP lower bound
        
#       standard rounding LambdaCC LP
        
        ub_lp_std_mean = round(Int64,mean(M["lp_std_ub_results"])) #standard rounding LambdaCC LP upper bound mean
        ub_lp_std_std = round(Int64,std(M["lp_std_ub_results"])) #standard rounding LambdaCC LP upper bound std
               
        rat_lp_std_mean = round(mean(M["lp_std_ub_results"] ./ M["lp_lb_results"]),sigdigits = dig2) #standard rounding LambdaCC LP ratio mean
        rat_lp_std_std = round(std(M["lp_std_ub_results"] ./ M["lp_lb_results"]),sigdigits = dig2) #standard rounding LambdaCC LP ratio std
       
        run_lp_std_mean = prun(mean(M["lp_std_tott_results"]),dig) #standard rounding LambdaCC total runtime mean
        run_lp_std_std = prun(std(M["lp_std_tott_results"]),dig) #standard rounding LambdaCC LP total runtime std
         
#         new rounding LambdaCC LP
        
        ub_lp_new_mean = round(Int64,mean(M["lp_new_ub_results"])) #new rounding LambdaCC LP upper bound mean
        ub_lp_new_std = round(Int64,std(M["lp_new_ub_results"]))#new rounding LambdaCC LP upper bound stdt
        
        rat_lp_new_mean = round(mean(M["lp_new_ub_results"] ./ M["lp_lb_results"]),sigdigits = dig2)#new rounding LambdaCC LP ratio mean
        rat_lp_new_std = round(std(M["lp_new_ub_results"] ./ M["lp_lb_results"]),sigdigits = dig2) #new rounding LambdaCC LP ratio std
        
        run_lp_new_mean = prun(mean(M["lp_new_tott_results"]),dig) #new rounding LambdaCC LP total runtime mean
        run_lp_new_std = prun(std(M["lp_new_tott_results"]),dig) #new rounding LambdaCC LP total runtime std
        
        
        status_lp = M["status_cclp"] # status--if false, we didn't converge 
#        
    catch
#       #Out of memory error"
        oom = "--"
        lb_cheaplp = oom
        ub_cheaplp_rand_std_mean = oom
        ub_cheaplp_rand_std_std = oom
        rat_cheaplp_rand_std_mean = oom
        rat_cheaplp_rand_std_std = oom
        run_cheaplp_rand_std_mean = oom
        run_cheaplp_rand_std_std = oom

        # lb_cheaplp = oom
        ub_cheaplp_rand_new_mean = oom
        ub_cheaplp_rand_new_std = oom
        rat_cheaplp_rand_new_mean = oom
        rat_cheaplp_rand_new_std = oom
        run_cheaplp_rand_new_mean = oom
        run_cheaplp_rand_new_std = oom


        status_cheaplp = 1.0
        
        lb_lp = oom
        ub_lp_std_mean= oom
        ub_lp_std_std = oom
        rat_lp_std_mean= oom
        rat_lp_std_std= oom
        run_lp_std_mean = oom
        run_lp_std_std = oom
        
        ub_lp_new_mean= oom
        ub_lp_new_std = oom
        rat_lp_new_mean= oom
        rat_lp_new_std= oom
        run_lp_new_mean = oom
        run_lp_new_std = oom
        
        status_lp = 1.0
    end

    oot = "Timeout"
   
    if status_cheaplp == 0.0

        lb_cheaplp = oot
        ub_cheaplp_rand_std_mean = oot
        ub_cheaplp_rand_std_std = oot
        rat_cheaplp_rand_std_mean = oot
        rat_cheaplp_rand_std_std = oot
        run_cheaplp_rand_std_mean = oot
        run_cheaplp_rand_std_std = oot

        lb_cheaplp = oot
        ub_cheaplp_rand_new_mean = oot
        ub_cheaplp_rand_new_std = oot
        rat_cheaplp_rand_new_mean = oot
        rat_cheaplp_rand_new_std = oot
        run_cheaplp_rand_new_mean = oot
        run_cheaplp_rand_new_std = oot
       
    end

    if status_lp == 0.0  

        lb_lp = oot
        ub_lp_std_mean= oot
        ub_lp_std_std = oot
        rat_lp_std_mean= oot
        rat_lp_std_std= oot
        run_lp_std_mean = oot
        run_lp_std_std = oot
        
        ub_lp_new_mean= oot
        ub_lp_new_std = oot
        rat_lp_new_mean= oot
        rat_lp_new_std= oot
        run_lp_new_mean = oot
        run_lp_new_std = oot
    end

    #--------------------------------------------------------------------------------------------------------------
    # print all values of cfp and lambdastc using standard and new rounding procedure
    # println("\\midrule")
    # if iscc_cheaplp == 1.0 && run_lp_std_mean != oot
    #     println("\\textsc{Lambda : $(lam)} & LB & $cfp_lb & $lb_cheaplp\$^*\$ & $lb_cheaplp\$^*\$ &$lb_lp & $lb_lp\\\\")
    # else
    #     println("\\textsc{Lambda : $(lam)} & LB & $cfp_lb & $lb_cheaplp  & $lb_cheaplp & $lb_lp & $lb_lp  \\\\")
    # end
    #     println("            & UB & $(ub_cfp_mean_round) \$\\tiny{\\pm $(ub_cfp_std_round)}\$ & $(ub_cheaplp_rand_std_mean) \$\\tiny{\\pm $(ub_cheaplp_rand_std_std)}\$ & $(ub_cheaplp_rand_new_mean) \$\\tiny{\\pm $(ub_cheaplp_rand_new_std)}\$ & $(ub_lp_std_mean) \$\\tiny{\\pm $(ub_lp_std_std)}\$ & $(ub_lp_new_mean) \$\\tiny{\\pm $(ub_lp_new_std)}\$  \\\\")
    #     println("\$n = $n \$ & Ratio & $rat_cfp_mean_round \$\\tiny{\\pm $rat_cfp_std_round}\$ & $rat_cheaplp_rand_std_mean \$\\tiny{\\pm $rat_cheaplp_rand_std_std}\$ & $rat_cheaplp_rand_new_mean \$\\tiny{\\pm $rat_cheaplp_rand_new_std}\$ & $rat_lp_std_mean \$\\tiny{\\pm $rat_lp_std_std}\$ &  $rat_lp_new_mean \$\\tiny{\\pm $rat_lp_new_std}\$ \\\\")
    #     println(" \$m = $m\$ & Run & $run_cfp_mean_round \$\\tiny{\\pm $run_cfp_std_round}\$ &  $run_cheaplp_rand_std_mean \$\\tiny{\\pm $run_cheaplp_rand_std_std}\$ &   $run_cheaplp_rand_new_mean \$\\tiny{\\pm $run_cheaplp_rand_new_std}\$ & $run_lp_std_mean \$\\tiny{\\pm $run_lp_std_std}\$ &  $run_lp_new_mean \$\\tiny{\\pm $run_lp_new_std}\$\\\\")
    
    #--------------------------------------------------------------------------------------------------------------
    #print only values of cfp and lambdastc using new rounding procedure
        # println("\\midrule")
    # if iscc_cheaplp == 1.0 && run_lp_std_mean != oot
    #     println("\\textsc{Lambda : $(lam)} & LB & $cfp_lb & $lb_cheaplp\$^*\$\\\\")
    # else
    #     println("\\textsc{Lambda : $(lam)} & LB & $cfp_lb & $lb_cheaplp \\\\")
    # end
    #     println("            & UB & $(ub_cfp_mean_round) \$\\tiny{\\pm $(ub_cfp_std_round)}\$ & $(ub_cheaplp_rand_mean) \$\\tiny{\\pm $(ub_cheaplp_rand_std)}\$   \\\\")
    #     println("\$n = $n \$ & Ratio & $rat_cfp_mean_round \$\\tiny{\\pm $rat_cfp_std_round}\$ & $rat_cheaplp_rand_mean \$\\tiny{\\pm $rat_cheaplp_rand_std}\$  \\\\")
    #     println(" \$m = $m\$ & Run & $run_cfp_mean_round \$\\tiny{\\pm $run_cfp_std_round}\$ &  $run_cheaplp_rand_mean \$\\tiny{\\pm $run_cheaplp_rand_std}\$ \\\\") 
   #---------------------------------------------------------------------------------------------------------------
   
   # Print table 

    println("\\midrule")
    if iscc_cheaplp == 1.0 && run_lp_std_mean != oot && (lam >= 0.5)
        println("            & $lam & $cfp_lb & $lb_cheaplp\$^*\$ & $(ub_cfp_mean_round) \$\\tiny{\\pm $(ub_cfp_std_round)}\$
        & $(ub_cheaplp_rand_std_mean) \$\\tiny{\\pm $(ub_cheaplp_rand_std_std)}\$ & $rat_cfp_mean_round \$\\tiny{\\pm $rat_cfp_std_round}\$
        & $rat_cheaplp_rand_std_mean \$\\tiny{\\pm $rat_cheaplp_rand_std_std}\$ & $run_cfp_mean_round \$\\tiny{\\pm $run_cfp_std_round}\$ &  $run_cheaplp_rand_std_mean \$\\tiny{\\pm $run_cheaplp_rand_std_std}\$ \\\\")
    else

        println("            & $lam & $cfp_lb & $lb_cheaplp & $(ub_cfp_mean_round) \$\\tiny{\\pm $(ub_cfp_std_round)}\$
        & $(ub_cheaplp_rand_new_mean) \$\\tiny{\\pm $(ub_cheaplp_rand_new_std)}\$ & $rat_cfp_mean_round \$\\tiny{\\pm $rat_cfp_std_round}\$
        & $rat_cheaplp_rand_new_mean \$\\tiny{\\pm $rat_cheaplp_rand_new_std}\$ & $run_cfp_mean_round \$\\tiny{\\pm $run_cfp_std_round}\$ &  $run_cheaplp_rand_new_mean \$\\tiny{\\pm $run_cheaplp_rand_new_std}\$ \\\\")

    end
        
end