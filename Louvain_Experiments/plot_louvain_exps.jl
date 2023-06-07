# using Plots
using MAT
using Statistics
using LaTeXStrings, Plots; 
gr()


L_m = Vector{Float64}()   # list of louvain ratio mean
L_std = Vector{Float64}()   # list of louvain ratio std
Lt_m = Vector{Float64}()  # list of louvain total runtimes mean
Lt_std = Vector{Float64}() # list of louvain total runtimes std

Lb = Vector{Float64}()  # list of lower bound
Lbt = Vector{Float64}()  #list of  lower bound time


CFP_m = Vector{Float64}()   # list of cfp ratio mean
CFP_std = Vector{Float64}()   #list of  cfp ratio std
CFPt_m= Vector{Float64}()  # list of cfp total runtimes mean
CFPt_std = Vector{Float64}()  #list of  cfp total runtimes std


S = Vector{Float64}()  #list of size of graph
V = Vector{Float64}()  #list of nunmber of nodes
E = Vector{Float64}()  #list of number of edges
Lambda = Vector{Float64}() #list of Lambda values


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
#Read SNAP graph
# graph = snaps[2]
# F = matread("../datafiles/simple-snap/simple-$graph")


#Read Facebook100 dataset
fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")
graph = fbs[8]
F = matread("../datafiles/Facebook100/$graph.mat")

A = Float64.(F["A"])
n = size(A,1)
lam_values = [0.4,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99] #

for lam in  lam_values

    F = matread("fb_louvain/$(graph)_louv_longround_lam_$(lam).mat")
#    F = matread("snap_louvain/$(graph)_louv_longround_lam_$(lam).mat")
    
   graph_stats = F["graph_stats"] # [n,m]
    
    lb = F["lb"]  #cfp lower bound
    lb_time = F["lb_time"]  #cfp lower bound time
 
    cl_ratio_mean = mean(F["cl_ratio_results"]) # mean of 15 cfp ratios
    cl_ratio_std = std(F["cl_ratio_results"])  # std of 15 cfp ratios
    
    cl_tott_mean = mean(F["cl_tott_results"]) # mean of 15 runtimes
    cl_tott_std = std(F["cl_tott_results"]) ## std of 15 runtimes
    
    llcl_ratio_mean = mean(F["llcl_ratio_results"])  # mean of 15 louvain ratios
    llcl_ratio_std = std(F["llcl_ratio_results"])  # std of 15 louvain ratios
    
    llcl_tott_mean = mean(F["llcl_tott_results"]) # mean of 15 louvain runtimes
    llcl_tott_std = std(F["llcl_tott_results"]) # std of 15 louvain runtimes
    
    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])
    push!(Lambda, graph_stats[3])

    push!(L_m,llcl_ratio_mean)
    push!(L_std,llcl_ratio_std)
    push!(Lt_m,llcl_tott_mean)
    push!(Lt_std,llcl_tott_std)
    
    push!(Lb,lb)
    push!(Lbt,lb_time)
   
    push!(CFP_m,cl_ratio_mean)
    push!(CFP_std,cl_ratio_std)
    push!(CFPt_m,cl_tott_mean)
    push!(CFPt_std,cl_tott_std)

end

## Plot approximations
lw = 1.5
ms = 6
gfs = 12
tfs = 10
titlesize = 14
s1 = 300
s2 = 200
sc = 2
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
color1 = :brown
color2 = :orange
color3 = :blue
ms1 = :star2
ms2 = :diamond
ms3 = :circle
title = ""
leg = :topright


lam_values_prun=[0.4,0.55,0.65,0.75,0.85,0.99]

# xlab = L"$\lambda$"
xlab = "Lambda"
ylab = "Approx Ratio"
ylims = (minimum(L_m)-0.2,maximum(CFP_m)+0.5)
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),color = :gray,xticks = (lam_values_prun,["$(v)" for v in lam_values_prun]))
title!(p,"cit-HepTh",titlefontsize = titlesize )
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs,ylim = ylims)

plot!(p,Lambda,[L_m L_m], fillrange=[L_m-L_std L_m+L_std], fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")
plot!(p,Lambda,[CFP_m CFP_m], fillrange=[CFP_m-CFP_std CFP_m+CFP_std], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")


plot!(p,Lambda,L_m, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "CFP+Louvain") #CFP+Louvain
plot!(p,Lambda,CFP_m, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "CFP") #CFP

savefig("Figures/fb_ratios_$(graph).pdf")
# savefig("Figures/snap_ratios_longround_$(graph).pdf")


##"0.4","0.55","0.65","0.75","0.85","0.95","0.99"
leg = :topright
ylab = "Runtime (sec)"
ylims = (minimum(CFPt_m)-maximum(CFPt_std)-0.3,maximum(Lt_m)+maximum(Lt_std)+1.0)
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :identity,legend = :topright,grid = false, size = (s1,s2),color = :gray,xticks = (lam_values_prun,["$(v)" for v in lam_values_prun]))
title!(p,"cit-HepTh",titlefontsize = titlesize)

xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs,ylim=ylims)

plot!(p,Lambda,[Lt_m Lt_m], fillrange=[Lt_m-Lt_std Lt_m+Lt_std], fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")
plot!(p,Lambda,[CFPt_m CFPt_m], fillrange=[CFPt_m-CFPt_std CFPt_m+CFPt_std], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")

plot!(p,Lambda,Lt_m, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "CFP+Louvain") #CFP+Louvain
plot!(p,Lambda,Lbt, fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "LB time") #LB time
plot!(p,Lambda,CFPt_m, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "CFP") #CFP

savefig("Figures/fb_runtimes_longround$(graph).pdf")
# savefig("Figures/snap_runtimes_longround_$(graph).pdf")
#-------------------------------------------------------------------------------------------------------------------------
