# FastLamCC
Faster approximation algorithms for LambdaCC

Code for faster approximation algorithms for LambdaCC by rounding the lower bounds into a clustering solution according. This code accompanies the anonymous paper submited at CIKM 2023.

We build on some of the previous code from the work of Veldt (ICML 2022). A copy of the code along with license is included in the include folder.

#Data
The datasets used are Facebook100 graphs and SNAP graphs that are publicly available and can be downloaded. SNAP graphs are standardized to make them undirected and remove self loops using 'standardize-snap-graphs.jl' 

#Experiments
There are three sets of experiments. 
1. CFP_LargeGraphs analyzes performance of the combinatorial CFP on large Facebook and SNAP graphs.

run_cfp_exps.jl calculates and stores results for each graph
plot_cfp_exps.jl plots the graphs based on the stored result files
Result files for Facebook and SNAP are stored in fb_cfp / snap_cfp respectively.

2. Louvain_Experiments includes results for a posteriori approximation guarantees using CFP

run_louvain_longround.jl and plot_louvain_exps.jl can be used reproduce the results for LambdaLouvain approximation with CFP lower bounds

fb_num_constraints/snap_num_constraints and plot_fb_constraints/plot_snap_constraints are used to plot the number of canonical LP, LambdaSTC LP and intermediate LP constraints for Facebook and SNAP graphs

'fb_louvain' and 'snap_louvain' folders have result files while the 'Figures' folder stores all the plots

3. LP_Experiments

run_lp_exps.jl file calculates LambdaCC clustering using CFP and cheaper LP method 
print_results.jl is used to generate the result values to create Table 2 in the arxiv paper

Folder 'src' includes the useful function definitions used in these experiments.
