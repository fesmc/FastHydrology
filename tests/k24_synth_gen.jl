# Generate the synthetic K24 test inputs once, so the Fortran driver and the
# Julia reference operate on bit-identical fields. Building them independently
# in each language makes them differ by ~1 ulp (Intel's libm vs Julia's for
# sin/cos/sqrt), which potential_filling's discrete local-minimum test can
# amplify into a large local difference -- an artifact of the test harness, not
# of either model.
using NCDatasets
const Nx, Ny = 48, 32
const dx, dy = 2000.0, 3000.0
A_const = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.0e-24

h=zeros(Nx,Ny); b=zeros(Nx,Ny); mask=zeros(Nx,Ny)
vb=zeros(Nx,Ny); A=zeros(Nx,Ny); mdot=zeros(Nx,Ny)
for j in 1:Ny, i in 1:Nx
    xr=(i-0.5*(Nx+1))/(0.5*Nx); yr=(j-0.5*(Ny+1))/(0.5*Ny); rr=sqrt(xr*xr+yr*yr)
    h[i,j]=2500.0*sqrt(max(0.0,1.0-rr))
    b[i,j]=-400.0+300.0*sin(3.0*pi*xr)*cos(2.0*pi*yr)
    mask[i,j]= h[i,j]>10.0 ? 1.0 : 0.0
    vb[i,j]=1.0e-6*(0.5+rr)
    A[i,j]=A_const*(1.0+0.5*cos(pi*xr))
    mdot[i,j]=1.0e-9*(1.0+0.5*sin(2.0*pi*xr)*sin(2.0*pi*yr))
end
NCDataset(ARGS[1],"c") do o
    defDim(o,"xc",Nx); defDim(o,"yc",Ny)
    defVar(o,"xc",collect(0:Nx-1).*dx,("xc",)); defVar(o,"yc",collect(0:Ny-1).*dy,("yc",))
    for (n,v) in [("h",h),("b",b),("mask",mask),("vb",vb),("A",A),("mdot",mdot)]
        defVar(o,n,v,("xc","yc"))
    end
end
println("wrote $(ARGS[1])  (A_const=$A_const)")
