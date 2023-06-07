using MAT
using LinearAlgebra
using SparseArrays
using Statistics
using Combinatorics
using Serialization
using Colors, ColorSchemes

include("../src/cfp_functions.jl")
include("../include/helpers.jl")
include("../include/faster_functions.jl")
using Plots

# Arrays to store the values
x_values = BigInt[]
num_owtri_values = BigInt[]
num_ow_values = BigInt[]
num_lp_values = BigFloat[]
class_values = Int[]
# fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")

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
    

data = Vector{Tuple{BigInt, BigInt, BigInt, BigFloat, Int}}()
for i in 1:33
    graph=snaps[i]
    loaded_dict = deserialize(open("snap_constraints/$(graph)_dictionary_data.jld"))

    push!(data, (loaded_dict["n"], loaded_dict["num_ow"],loaded_dict["num_tri"], loaded_dict["lp_constr"], loaded_dict["classnum"]))
end


sort!(data, by = x -> x[1])
x_values = [x[1] for x in data]
num_ow_values = [x[2] for x in data]
num_owtri_values = [(x[2]+x[3]) for x in data]
num_lp_values = [x[4] for x in data]
class_values = [x[5] for x in data]

# # Plotting
lw = 1.5
ms = 6
gfs = 12
tfs = 10
titlesize = 14
s1 = 300
s2 = 225 #225
sc = 2
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
color1 = :brown
color2 = :green
color3 = :purple
ms1 = :star2
ms2 = :diamond
ms3 = :circle
title = ""
leg = :topright

xlab = "Number of nodes"
ylab = "Constraint size"

# Define a color palette with different colors for each class

color_palette = [RGB(0.2,0.133,0.533),RGB(0.533,0.8,0.933),RGB(0.267,0.667,0.6),RGB(0.067,0.467,0.2),RGB(0.6,0.6,0.2),RGB(0.867,0.8,0.467),RGB(0.8,0.4,0.467),RGB(0.533,0.133,0.333)]

label_values = ["loc-social","o-social","web","comm","road","prod","collab","cit"]

# Plotting
f = plot(
    title = "",
    fontfamily = "helvetica",
    linewidth = lw,
    size = (s1, s2),
    titlefontsize = titlesize,
    grid = false,
    yscale = :log10,
    xscale = :log10,
    legend = :topright,
    legendfontsize = 6,
    legend_background_color = colorant"rgba(255,255,255,0.6)"

)

plot!(f, x_values, num_lp_values, line = :solid, linewidth = lw, label = "", color = :black) # canonical LP


xaxis!(f, xlab, tickfontsize = tfs, guidefontsize = gfs)
yaxis!(f, ylab, tickfontsize = tfs, guidefontsize = gfs)

# Update the legend with class labels
for class_val in 1:8
    label="--"
    if(class_val == 1)
        label = "loc-social"   
    elseif(class_val == 2)
        label = "o-social" 
    elseif(class_val == 3)
        label = "web" 
    elseif(class_val == 4)
        label = "comm" 
    elseif(class_val == 5)
        label = "road" 
    elseif(class_val == 6)
        label = "prod" 
    elseif(class_val == 7)
        label = "collab" 
    elseif(class_val == 8)
        label = "cit" 
    else
        println("label not found")
    end
 
    color = color_palette[class_val]


    S = findall(x->x==class_val,class_values)

    scatter!(f,x_values[S],num_owtri_values[S],markershape=ms1 , markerstrokewidth = 0, label = "", color = color_palette[class_val]) #LambdaSTC LP

    scatter!(f,x_values[S],num_ow_values[S],markerstrokewidth = 0, label = label, color = color_palette[class_val]) #LambdaSTC LP

    # scatter!(f, [NaN], [NaN], markerstrokewidth = 0, markercolor = color_palette[class_val], label = label)
end
ylims!(f, ylims(f)[1], ylims(f)[2]*100)
# xlims!(f, xlims(f)[1], xlims(f)[2]*10)


# savefig("Figures/facebook_num_constraints.pdf")
savefig("Figures/snap_num_all_constraints.pdf")
