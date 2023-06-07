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
num_tri_values = BigInt[]
num_lp_values = BigFloat[]
fbs = readlines("../datafiles/fb-graph-names2-Copy1.txt")

data = Vector{Tuple{BigInt, BigInt, BigInt, BigFloat, Int}}()

for i in 1:100
    graph=fbs[i]
    loaded_dict = deserialize(open("fb_constraints/$(graph)_dictionary_data.jld"))
    push!(data, (loaded_dict["n"], loaded_dict["num_ow"],loaded_dict["num_tri"], loaded_dict["lp_constr"], 1))
end

sort!(data, by = x -> x[1])
x_values = [x[1] for x in data]
num_ow_values = [x[2] for x in data]
num_tri_values = [x[3] for x in data]
num_owtri_values = [(x[2]+x[3]) for x in data]
num_lp_values = [x[4] for x in data]
class_values = [x[5] for x in data]


# Plotting
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

#Plotting graphs-------------------------------------------------------------------------------------
f = plot(title = "",fontfamily = "helvetica",linewidth = lw, size = (s1,s2), titlefontsize = titlesize, grid = false, yscale = :log10, xscale= :log10 ,legend = :topleft, legend_background_color = colorant"rgba(255,255,255,0.6)")
scatter!(f,x_values,num_owtri_values,markershape = ms1, markerstrokewidth = 0, label = "intermed. LP", color = :maroon) #intermed LP
scatter!(f,x_values,num_ow_values,markerstrokewidth = 0, label = "LambdaSTC LP", color = :green) #LambdaSTC LP
plot!(f, x_values, num_lp_values, line = :solid, linewidth = lw, label = "canonical LP", color = :black) # canonical LP

xaxis!(f,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(f,ylab,tickfontsize=tfs,guidefontsize = gfs)
# ylims!(f, ylims(f)[1], ylims(f)[2]*10)
savefig("Figures/facebook_num_all_constraints.pdf")
# savefig("Figures/snap_num_constraints.pdf")
