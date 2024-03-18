include("../Functions/Dependencies_3d.jl");
using MKLSparse

β = 100;
Ω = 0.
N  = 8;
ℓ = 1;
box_size = 6; #will give a cube domain from -box to box
ϵ = 1/2;
mesh = Mesh3d.generate(N,box_size,ℓ)
sparsity = Assemble.MatrixSparsity(mesh);

tid_all = time();


@time ∇vᵢ∇vⱼ  = Assemble.∇vᵢ∇vⱼ(mesh,Quad4,sparsity)
v(x) = 1.;
@time vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad6,sparsity);

V(x) = 1*dot(x,x)+100*dot(sin.(pi*x),sin.(pi*x));

dynamics = true
Vₜ(x) = 2*dot(x,x);



tid = time();
Vvᵢvⱼ = Assemble.vᵢvⱼ(mesh,V,Quad8,sparsity);
println("computed potential ", time()-tid);
#vᵢLzvⱼ = Assemble.vᵢLzvⱼ(mesh,Quad6);
P = Assemble.Coarse2Fine(mesh);



#B = ∇vᵢ∇vⱼ ;

tid = time()
ϕ = Assemble.ϕ_canonical(mesh,∇vᵢ∇vⱼ,vᵢvⱼ,P);
println("computed SLOD-basis ", time()-tid);

#b = load("Matrices.jld");
#ϕ = b["ϕ"];
#Vvᵢvⱼ = b["Vvᵢvⱼ"]#;,ϕ, "Vvᵢvⱼ",Vvᵢvⱼ);

#save("Matrices.jld","ϕ",ϕ,"Vvᵢvⱼ",Vvᵢvⱼ)
b = 0;
ϕ = ϕ[:,mesh.dofs];

∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ[mesh.dofs,mesh.dofs]*ϕ')
φᵢφⱼ = ϕ*(vᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');  φᵢφⱼ +=φᵢφⱼ'; φᵢφⱼ /=2;
Vφᵢφⱼ  = ϕ*(Vvᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');
#φᵢLzφⱼ  = ϕ*(vᵢLzvⱼ[mesh.dofs,mesh.dofs]*ϕ');

∇vᵢ∇vⱼ = 0;
Vvᵢvⱼ = 0; 

SpaceDim_C = size(ϕ,1);

tid = time();
ωₕ=ω_module.PreAllocateW(vᵢvⱼ[mesh.dofs,mesh.dofs])
ω_module.Compute_Wh(mesh,ωₕ,Quad8);

ω̃ₕ  = ω_module.Compute_tilde(ωₕ )

#---------------------------------Initial guess -----------------------------------------------#


UGS = Assemble.initial_guess_IterativeSolvers(Ω,mesh,φᵢφⱼ,vᵢvⱼ,mesh.dofs,ϕ)

println("time all pre comp ", time()-tid_all);
#----------------------------------------------------------------------------------------------#
max_it = 1000
tol_res = 0.05 #shift to inverse iteration
tol_stop =   1e-6;
UGS, E,conv = J_METHOD_ωₕ_IterativeSolvers(mesh,vᵢvⱼ[mesh.dofs,mesh.dofs],∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,UGS,ϕ,max_it,ωₕ,ω̃ₕ ,ϵ,β,Quad8,tol_res,tol_stop)	
save("./Results/J_method_smooth_N_"*string(N)*"_beta_"*string(β)*"_3D_new_basis.jld","UGS",UGS,"E",E,"conv",conv,"t",time()-tid)
ω̃ₕ = 0;
tid = time()
if(dynamics)
	Vₜvᵢvⱼ = Assemble.vᵢvⱼ(mesh,Vₜ,Quad8,sparsity);
	Vₜφᵢφⱼ = ϕ*(Vₜvᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');
	println("Computed smooth part in ", time()-tid);


	T = 1; Nt = 2^7; max_it = 20; save_it = 1; save_here = "./Results/"; TOL = 10^-8;
	CG_q2_ωₕ_IterativeSolvers(UGS,T,Nt,ϕ,∇φᵢ∇φⱼ,φᵢφⱼ,Vₜφᵢφⱼ,ωₕ,TOL,max_it,β,ϵ,save_it,save_here,Quad8)
end

#---------------------------------------------------------------------------------------------#
#=
using PlotlyJS
UGS_h = ϕ'*UGS;

UGS_plot = 0im*zeros(size(mesh.p,2));
UGS_plot[mesh.dofs] .= UGS_h;
UGS_plot = reshape(UGS_plot,(3*N+1),(3*N+1),(3*N+1))
data = range(-mesh.box_size, stop=mesh.box_size, length=3*N+1)

X, Y, Z = mgrid(data, data, data)
ρ = abs.(UGS_plot).^2
#ρ = sin.(X .* Y .* Z) ./ (X .* Y .* Z)


plot(volume(

    x=X[:],

    y=Y[:],

    z=Z[:],

    value= ρ[:],

    isomin=0.01,

 #   isomax=0.8,

    opacity=0.1, # needs to be small to see through all surfaces

    surface_count=17, # needs to be a large number for good volume rendering

))



=#

