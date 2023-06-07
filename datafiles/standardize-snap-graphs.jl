using MAT
using LinearAlgebra
using SparseArrays
# Download and save all of the following graphs in .mat form in the snap graphs folder.
# .mat files for these graphs can all be found on the SuiteSparse Matrix Collection 
# https://sparse.tamu.edu/
snaps = [
    "amazon0302.mat";
    "amazon0312.mat";
    "amazon0505.mat";
    "amazon0601.mat";
    "ca-AstroPh.mat" ;
    "ca-CondMat.mat";
    "ca-GrQc.mat";
    "ca-HepPh.mat" ;
    "ca-HepTh.mat";
    "cit-HepPh.mat";
    "cit-HepTh.mat" ;
    "cit-Patents.mat";
    "com-Amazon.mat";
    "com-DBLP.mat";
    "com-LiveJournal.mat";
    "com-Youtube.mat";
    "email-Enron.mat";
    "email-EuAll.mat";
    "loc-Brightkite.mat";
    "loc-Gowalla.mat";
    "roadNet-CA.mat";
    "roadNet-PA.mat";
    "roadNet-TX.mat";
    "soc-Epinions1.mat";
    "soc-LiveJournal1.mat";
    "soc-Slashdot0811.mat";
    "soc-Slashdot0902.mat";
    "web-BerkStan.mat";
    "web-Google.mat";
    "web-NotreDame.mat";
    "web-Stanford.mat";
    "wiki-Talk.mat";
    "wiki-topcats.mat"
]


# Update this path to wherever you store the snap networks
pathtosnaps = "../snap/"

# This code standardizes the graphs so that they are undirected, unweighted, and have no self-loops

for graph = snaps

    F = matread("snap/$(graph)")
    P = F["Problem"]
    A = P["A"]

    n = size(A,1)

    # Remove self-loops
    for i = 1:n
        A[i,i] = 0  
    end
    dropzeros!(A)

    println("$graph $n")

    # Remove edge directions
    A = sparse(A+A')

    # Remove weights
    I,J,V = findnz(A)
    A = sparse(I,J,1.0,n,n)

    matwrite("simple-snap/simple-$graph", Dict("A"=>A);compress=true)

    @assert(SparseArrays.issymmetric(A))
end
