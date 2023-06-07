using Random
using SparseArrays
using LinearAlgebra
using MAT

function coverLabel(A,lam)
    
    gam=lam/(1-lam)
    
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    Add = Dict()
    Del = Dict()
    Res = Dict()  # or can use adjacency matrix as Residual matrix
    sumY = 0
    for i = 1:n
        N = Neighbs[i]
  
        for jj = 1:length(N)
            j = N[jj]
            ij1 = max(i,j)
            ij2 = min(i,j)
            if haskey(Res,(ij1,ij2)) && Res[(ij1,ij2)]==0
                continue
            end

            for kk = jj+1:length(N)
                k = N[kk]
                ik1 = max(i,k)
                ik2 = min(i,k)

                jk1 = max(j,k)
                jk2 = min(j,k)

                if A[k,j] > 0
                    continue
                end
#             if (haskey(Res,(ij1,ij2)) && Res[(ij1,ij2)]==0) 
#                 continue
#             end
            
            
                if !haskey(Res,(ij1,ij2))
                    Res[(ij1,ij2)]=(1-lam) # or just 1
                end 
                if !haskey(Res,(ik1,ik2))
                    Res[(ik1,ik2)]=(1-lam) # or just 1
                end
                if !haskey(Res,(jk1,jk2))
                    Res[(jk1,jk2)]=lam # or just 1
                end
            
                if Res[(ij1,ij2)]>0 && Res[(ik1,ik2)]>0 && Res[(jk1,jk2)]>0
                    y = min(Res[(ij1,ij2)],Res[(ik1,ik2)],Res[(jk1,jk2)])
                    Res[(ij1,ij2)] = Res[(ij1,ij2)] - y
                    Res[(ik1,ik2)] = Res[(ik1,ik2)] - y
                    Res[(jk1,jk2)] = Res[(jk1,jk2)] - y
                    sumY +=y
                
                
                end
            
                if Res[(ij1,ij2)] ==0
                    break
                end
                
            
                    # if we reach here: (i,j) and (i,k) are edges, but (j,k) is not
                    # So (i,j,k) is an open wedge centered at i.
                    # Also, neither (i,j) nor (i,k) nor (j,k) already appears in our open wedge set.

#                     Del[(ik1,ik2)] = true
#                     Del[(ij1,ij2)] = true
#                     Add[(jk1,jk2)] = true

                    # Now exit this and start with a new j, as we have already
                    # used edge (i,j) and we don't need to check for more
                    # edges involving (i,j)
#                 break
            
            end
        end
    end
    
    for ((k1,k2),v) in Res
        if(v==0 && A[k1,k2]==1 )
            Del[(k1,k2)] = true
        end
        if (v==0 && A[k1,k2]==0)
            Add[(k1,k2)] = true
        end
        
    end
    
    Edel = collect(keys(Del))
    Eadd = collect(keys(Add))
        
#     bd = ((1-lam)*(size(Edel,1))+((lam)*size(Eadd,1)))   
    
    return Edel, Eadd,round(Int64,sumY) # round(Int64,bd)

    
end