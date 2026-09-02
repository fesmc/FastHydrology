using NCDatasets, Statistics, Printf
f=NCDataset(ARGS[1]); j=NCDataset(ARGS[2])
vars = length(ARGS)>2 ? split(ARGS[3],",") : ["phi0","absg","gsx","gsy","absgs","q","W","N","p_w"]
m = Bool.(j["mask"][:,:] .== 1.0)
@printf("%-8s %14s %14s %12s   %s\n","field","max|J|","max|F-J|","rel_max","(all cells rel_max)")
for v in vars
    F=Float64.(f[v][:,:]); J=Float64.(j[v][:,:])
    Fm,Jm = F[m],J[m]
    s=maximum(abs,Jm); d=maximum(abs,Fm.-Jm)
    sa=maximum(abs,J); da=maximum(abs,F.-J)
    @printf("%-8s %14.6e %14.6e %12.3e   %12.3e\n",v,s,d,s>0 ? d/s : 0.0, sa>0 ? da/sa : 0.0)
end
close(f);close(j)
