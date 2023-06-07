# Plotting approximation ratios and runtimes for three graphs
using MAT
using Statistics
using LaTeXStrings, Plots; 
gr()

# Result for graph1
Lb = Vector{Float64}()  # list of lower bounds
Lbt = Vector{Float64}()  # list of lower bound time


CFP_m = Vector{Float64}()   # list of cfp ratio mean
CFP_std = Vector{Float64}()   # list of cfp ratio std
CFPt_m= Vector{Float64}()  # list of cfp total runtimes mean
CFPt_std = Vector{Float64}()  # list of cfp total runtimes std


S = Vector{Float64}() #list of graph size (S=V+E)
V = Vector{Float64}() #list of node size
E = Vector{Float64}() #list of edge size
# =---------------------------------------------------------
#Results for graph 2
Lb2 = Vector{Float64}()  # lower bound
Lbt2 = Vector{Float64}()  # lower bound time


CFP_m2 = Vector{Float64}()   # cfp ratio mean
CFP_std2 = Vector{Float64}()   # cfp ratio std
CFPt_m2= Vector{Float64}()  # mfp total runtimes mean
CFPt_std2 = Vector{Float64}()  # mfp total runtimes std


S2 = Vector{Float64}()
V2 = Vector{Float64}()
E2 = Vector{Float64}()
# =---------------------------------------------------------
#Results for graph 3
Lb3 = Vector{Float64}()  # lower bound
Lbt3 = Vector{Float64}()  # lower bound time


CFP_m3 = Vector{Float64}()   # cfp ratio mean
CFP_std3 = Vector{Float64}()   # cfp ratio std
CFPt_m3= Vector{Float64}()  # mfp total runtimes mean
CFPt_std3 = Vector{Float64}()  # mfp total runtimes std


S3 = Vector{Float64}()
V3 = Vector{Float64}()
E3 = Vector{Float64}()
# ----------------------------------------------------------
Lambda = Vector{Float64}()

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

fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")


graph = fbs[64]
# graph = snaps[5]

graph2 = snaps[12]
# graph2 = fbs[48]

graph3 = snaps[21]
# graph = fbs[3]

F = matread("../datafiles/Facebook100/$graph.mat")
# F = matread("../datafiles/simple-snap/simple-$graph.mat")

F2 = matread("../datafiles/simple-snap/simple-$graph2.mat")
# F2 = matread("../datafiles/Facebook100/$graph2.mat") 

F3 = matread("../datafiles/simple-snap/simple-$graph3.mat")
# F3 = matread("../datafiles/Facebook100/$graph3.mat") 

lam_values = [0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99]
for lam in  lam_values

   #get results for first graph
    # F = matread("snap_cfp/$(graph)_snap_cfp_lam_$(lam).mat")
    F = matread("fb_cfp/$(graph)_fb_cfp_lam_$(lam).mat")

    graph_stats = F["graph_stats"] # [n,m,lam]
    
    lb = F["lb"] #lower bound
    lb_time = F["lb_time"] #lower bound time
 
    cl_ratio_mean = mean(F["cl_ratio_results"]) #cfp ratio mean
    cl_ratio_std = std(F["cl_ratio_results"])#cfp ratio std
    
    cl_tott_mean = mean(F["cl_tott_results"]) #cfp totaltime mean
    cl_tott_std = std(F["cl_tott_results"]) #cfp totaltime std
     
    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])
    push!(Lambda, graph_stats[3])
 
    push!(Lb,lb)
    push!(Lbt,lb_time)
   
    push!(CFP_m,cl_ratio_mean)
    push!(CFP_std,cl_ratio_std)
    push!(CFPt_m,cl_tott_mean)
    push!(CFPt_std,cl_tott_std)

    # get results for second graph--------------------------------------------------
    # F2 = matread("fb_cfp/$(graph2)_fb_cfp_lam_$(lam).mat")
    F2 = matread("snap_cfp/$(graph2).mat_snap_cfp_lam_$(lam).mat")

    graph_stats2 = F2["graph_stats"] # [n,m,open wedges, triangles]
    
    lb2 = F2["lb"]
    lb_time2 = F2["lb_time"]
 
    cl_ratio_mean2 = mean(F2["cl_ratio_results"])
    cl_ratio_std2 = std(F2["cl_ratio_results"])
    
    cl_tott_mean2 = mean(F2["cl_tott_results"])
    cl_tott_std2 = std(F2["cl_tott_results"])
    
    push!(V2, graph_stats2[1])
    push!(E2, graph_stats2[2])
    push!(S2, graph_stats2[1]+graph_stats2[2])
    
    push!(Lb2,lb2)
    push!(Lbt2,lb_time2)
   
    push!(CFP_m2,cl_ratio_mean2)
    push!(CFP_std2,cl_ratio_std2)
    push!(CFPt_m2,cl_tott_mean2)
    push!(CFPt_std2,cl_tott_std2)

    # get results for  third graph--------------------------------------------------
    # F2 = matread("fb_cfp/$(graph2)_fb_cfp_lam_$(lam).mat")
    F3 = matread("snap_cfp/$(graph3)_snap_cfp_lam_$(lam).mat")
    graph_stats3 = F3["graph_stats"] # [n,m,lam]
    
    lb3 = F3["lb"]
    lb_time3 = F3["lb_time"]
 
    cl_ratio_mean3 = mean(F3["cl_ratio_results"])
    cl_ratio_std3 = std(F3["cl_ratio_results"])
    
    cl_tott_mean3 = mean(F3["cl_tott_results"])
    cl_tott_std3 = std(F3["cl_tott_results"])
  
    
    push!(V3, graph_stats3[1])
    push!(E3, graph_stats3[2])
    push!(S3, graph_stats3[1]+graph_stats3[2])
    
    push!(Lb3,lb3)
    push!(Lbt3,lb_time3)
   
    push!(CFP_m3,cl_ratio_mean3)
    push!(CFP_std3,cl_ratio_std3)
    push!(CFPt_m3,cl_tott_mean3)
    push!(CFPt_std3,cl_tott_std3)

end

## Plot approximations
lw = 1.5
ms = 6
gfs = 12 #12
tfs = 10
titlesize = 14
s1 = 300
s2 = 200
sc = 2
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
color1 = :indigo
color2 = :sienna1
color3 = :cyan4
ms1 = :star2
ms2 = :diamond
ms3 = :circle
title = ""
leg = :topright

lam_values_prun=[0.55,0.65,0.75,0.85,0.9,0.99]

#Plot approximation ratios
# xlab = L"$\lambda$"
xlab = "Lambda"
ylab = "Approx Ratio"
ylims = (min(minimum(CFP_m),minimum(CFP_m2),minimum(CFP_m3))-0.2,max(maximum(CFP_m),maximum(CFP_m2),maximum(CFP_m3))+1.0)
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),legendfontsize = 7,color = :gray,xticks = (lam_values_prun,["$(v)" for v in lam_values_prun]))
title!(p,"CFP Ratios",titlefontsize = titlesize) #,ylim = ylims)
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)


plot!(p,Lambda,[CFP_m CFP_m], fillrange=[CFP_m-CFP_std CFP_m+CFP_std], fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")

plot!(p,Lambda,CFP_m, fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "Texas84") 

# Plot graph 2-------------------------------------------------------------------------------
plot!(p,Lambda,[CFP_m2 CFP_m2], fillrange=[CFP_m2-CFP_std2 CFP_m2+CFP_std2], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")

plot!(p,Lambda,CFP_m2, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "cit-Patents")
# Plot graph 3-------------------------------------------------------------------------------
plot!(p,Lambda,[CFP_m3 CFP_m3], fillrange=[CFP_m3-CFP_std3 CFP_m3+CFP_std3], fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "")

plot!(p,Lambda,CFP_m3, fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "roadnet-CA")

# savefig("Figures/snap_ratios_$(graph).pdf")
savefig("Figures/fbsnap_ratios_$(graph)_$(graph2)_$(graph3).pdf")

# RUNTIME===================================================================================
##"0.4","0.55","0.65","0.75","0.85","0.95","0.99"
leg = :topright
ylab = "Runtime (sec)"
ylims = (min(minimum(CFPt_m),minimum(CFPt_m2))-0.2,max(maximum(CFPt_m),maximum(CFPt_m2))+1.0)
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :log,legend = leg,grid = false, size = (s1,s2),color = :gray,xticks = (lam_values_prun,["$(v)" for v in lam_values_prun]))
title!(p,"CFP Runtimes ",titlefontsize = titlesize)

xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)


plot!(p,Lambda,[CFPt_m CFPt_m], fillrange=[CFPt_m-CFPt_std CFPt_m+CFPt_std], fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")

plot!(p,Lambda,CFPt_m, fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")

# Plot graph 2 runtimes=======================================================================================================================================
plot!(p,Lambda,[CFPt_m2 CFPt_m2], fillrange=[CFPt_m2-CFPt_std2 CFPt_m2+CFPt_std2], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")

plot!(p,Lambda,CFPt_m2, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")

# Plot graph3 runtimes=======================================================================================================================================
plot!(p,Lambda,[CFPt_m3 CFPt_m3], fillrange=[CFPt_m3-CFPt_std3 CFPt_m3+CFPt_std3], fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
markershape = ms3, markerstrokecolor = color3, label = "")

plot!(p,Lambda,CFPt_m3, fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
markershape = ms3, markerstrokecolor = color3, label = "")

# savefig("Figures/fb_runtimes_$(graph).pdf")
savefig("Figures/fbsnap_runtimes_$(graph)_$(graph2)_$(graph3).pdf")
#-------------------------------------------------------------------------------------------------------------------------