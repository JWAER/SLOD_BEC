include("../Functions/Dependencies_3d.jl");


prealloc = 2*10^8;
#----------------------------Generate or load square domain -------------------------------------#
ℓ = 1;
Nh = 1; #refinement, for smooth problems Nh = 1;
N = 30;
β = 50
Ω = 0.0;
ϵ = 1/2; 
box = 3; #will give a square domain from -box to box

mesh = Mesh3d.generate(N,Nh,box,ℓ);
sparsity = Assemble.MatrixSparsity(mesh);

Vd(x) =2*mod(floor(x[1])+floor(x[2])+floor(x[3]),2);
 #therefore only compute if |xᵢ|<ℓH

γ_x = 1.0;
γ_y = 1.0;
Vs(x) = 1/2*dot(x,x);

dynamics = true;
Vₜ(x) = 1/2*dot(x,x);
#------------------------------Compute All Matrix Quantities------------------------------------#
tid_pre = time();
	sparsity = Assemble.MatrixSparsity(mesh)

	tid = time();
	∇vᵢ∇vⱼ = Assemble.∇vᵢ∇vⱼ(mesh,Quad4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	v(x) = 1.0
	vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	tid=time();
	#vᵢvⱼVd = Assemble.vᵢvⱼ(mesh,Vd,Quad_9,sparsity);
	vᵢvⱼVd = Assemble.vᵢvⱼ(mesh,Vd,Quad6,sparsity); #piecewise constant
	println("Computed rough matrix in ", time()-tid);

	tid = time();
	vᵢvⱼVs = Assemble.vᵢvⱼ(mesh,Vs,Quad8,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.Assemble_P(mesh)
	println("Computed P in ", time()-tid);

tid_pre = time()-tid_pre;



vᵢvⱼV = vᵢvⱼVs+vᵢvⱼVd;
Aᵥ = ϵ*∇vᵢ∇vⱼ+vᵢvⱼVd;

tid = time();
ϕ = Assemble.ϕ_refinement2(mesh,Aᵥ,vᵢvⱼ,P,ℓ,sparsity)#,∇vᵢ∇vⱼ); #assembles EVERY SLOD-function ... possibly ok in 2d
tid_SLOD = time()-tid;

#------------------------- Restrict to H^1_0   ----------------------------------
dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);

ϕ = sparse(ϕ)[:,dofs_f];
∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ[dofs_f,dofs_f]*ϕ')
φᵢφⱼ = ϕ*(vᵢvⱼ[dofs_f,dofs_f]*ϕ')
Vφᵢφⱼ = ϕ*(vᵢvⱼV[dofs_f,dofs_f]*ϕ');

SpaceDim_C = size(ϕ,1);
tid = time();
ωₕ=ω_module.PreAllocateW(vᵢvⱼ[mesh.dofs,mesh.dofs])
ω_module.Compute_Wh(mesh,ωₕ,Quad8);
ω̃ₕ  = ω_module.Compute_tilde(ωₕ,length(ωₕ.Val)   )


#-----------------------------------Solve for minimizier --------------------------------------#
UGS = Assemble.initial_guess_IterativeSolvers(Ω,mesh,φᵢφⱼ,vᵢvⱼ,mesh.dofs,ϕ)

max_it = 100
tol_res = 0.1
UGS, E,conv = J_METHOD_ωₕ_IterativeSolvers(mesh,vᵢvⱼ[mesh.dofs,mesh.dofs],∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,UGS,ϕ,max_it,ωₕ,ω̃ₕ ,ϵ,β,Quad8,0.05)	

save("../Results/J_method_smooth_N_"*string(N)*"_beta_"*string(β)*"_3D_new_basis.jld","UGS",UGS,"E",E,"conv",conv,"t",tid,"t_pre",tid_pre,"t_SLOD",tid_SLOD)

tid = time()
if(dynamics)
Vₜvᵢvⱼ = Assemble.vᵢvⱼ(mesh,Vₜ,Quad8,sparsity);
	Vₜφᵢφⱼ = ϕ*(Vₜvᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');
	println("Computed smooth part in ", time()-tid);


	T = 10; Nt = 2^7; max_it = 20; save_it = 1; save_here = "./Results/Final/"; TOL = 10^-8;
	CG_q3_ωₕ_IterativeSolvers(UGS,T,Nt,ϕ,∇φᵢ∇φⱼ,φᵢφⱼ,Vₜφᵢφⱼ,ωₕ,TOL,max_it,β,ϵ,save_it,save_here)

end





#end

#-----------------PLOT on cartesian grid ---------------------------------------#
#=
using PlotlyJS
ugs = 0im*zeros(size(mesh.p,2));
ugs[dofs_f] = ϕ'*UGS;
ugs = reshape(ugs,3N*Nh+1,3N*Nh+1);


plot(contour(

   # colorscale="Jet",
    x=mesh.p[1,1:3N*Nh+1],
    y=mesh.p[1,1:3N*Nh+1], 

    z=abs.(ugs).^2

))

#plot(surface(z=ugs.^2))
=#
#=
plot(heatmap(z=abs.(ugs).^2,
    x=mesh.p[1,1:3N*Nh+1],
    y=mesh.p[1,1:3N*Nh+1]))


=#





