include("../Functions/Dependencies_2d.jl");


prealloc = 2*10^8;
#----------------------------Generate or load square domain -------------------------------------#
ℓ = 2;
for Nh = [2]; #refinement, for smooth problems Nh = 1;
Ns = [12 18 24 36 48 72 96 144 192 384]

#for it = 3
#for it = 2
for it = 9
βs = 100;
Ω = 0.0;
ϵ = 1/2; 
box = 12; #will give a square domain from -6 to 6
N = Ns[it]
Vd_refinement_factor = 1#Vd_refinement_factors[it];

mesh = SLOD_MESH.Generate(N,Nh,box,ℓ);

Vd(x) = 2*sum(x.>0);#5 + min.(floor.(2*sin(π*x[1]/3)*sin(π*x[2]/3)),1);
#Vd(x) = (x[1]>0);
#Vd(x) = 2*abs(x[1]/epsilon-round(x[1]/epsilon))+2*abs(x[2]/epsilon-round(x[2]/epsilon))
γ_x = 1.0;
γ_y = 1.0;
Vs(x) = 1/2*dot(x,x);

dynamics = false;
βt = 20;
Vₜ(x) = 5 + min.(floor.(2*sin(π*x[1]/1)*sin(π*x[2]/1)),1);
#------------------------------Compute All Matrix Quantities------------------------------------#
tid = time();
	sparsity = Assemble.MatrixSparsity(mesh);

	tid = time();
	∇vᵢ∇vⱼ = Assemble.∇vᵢ∇vⱼ(mesh,Quad_4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad_6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	vᵢvⱼVd = Assemble.Vvᵢvⱼ(Vd,mesh,Quad_6,sparsity);

	tid = time();
	vᵢvⱼVs = Assemble.Vvᵢvⱼ(Vs,mesh,Quad_9,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.P(mesh)
	println("Computed P in ", time()-tid);



vᵢvⱼV = vᵢvⱼVs+vᵢvⱼVd;
Aᵥ = ϵ*∇vᵢ∇vⱼ+vᵢvⱼVd;

tid = time();
ϕ = Assemble.ϕ(mesh,Aᵥ,vᵢvⱼ,P,ℓ); #

tid_assembly = time()-tid;
#ϕ = Assemble.ϕ2(mesh,Aᵥ,vᵢvⱼ,P,ℓ); #s

println("Assemby  ", tid_assembly);

#------------------------- Restrict to H^1_0   ----------------------------------
	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);

	ϕ = sparse(ϕ)[:,dofs_f];

	∇vᵢ∇vⱼ = ∇vᵢ∇vⱼ[dofs_f,dofs_f]
	vᵢvⱼ = vᵢvⱼ[dofs_f,dofs_f];
	vᵢvⱼV =vᵢvⱼV[dofs_f,dofs_f]
	vᵢLzvⱼ =vᵢLzvⱼ[dofs_f,dofs_f];
	

	∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ*ϕ')
	φᵢφⱼ = ϕ*(vᵢvⱼ*ϕ')
	Vφᵢφⱼ = ϕ*(vᵢvⱼV*ϕ');
	φᵢLzφⱼ = ϕ*(vᵢLzvⱼ*ϕ');
	

#------------------------------------------------------------------------------------------------#

#---------------------------------------- Compute Tensor ----------------------------------------#

	tid = time();
	ωₕ=ω_module.PreAllocateW(vᵢvⱼ[dofs_f,dofs_f])
	ω_module.Compute_Wh(mesh,ωₕ);

	ω̃ₕ  = ω_module.Compute_tilde(ωₕ)

#----------------------------------------------------------------------------------------------#

#-----------------------------------Solve for minimizier --------------------------------------#


φᵢφⱼ_lu = lu(φᵢφⱼ); # compute lu once
UGS = Assemble.initial_guess(Ω,mesh,φᵢφⱼ_lu,vᵢvⱼ,dofs_f,ϕ)

#----------------------------------------------------------------------------------------------#
#max_it = 20
tid = time()
max_it = 1000
tol_Res = .1;
tol_stop = 1e-7
UGS, Ẽ,E,conv = J_METHOD_ωₕ(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,UGS,ϕ,max_it,φᵢφⱼ_lu,ωₕ,ω̃ₕ ,ϵ,βs,Quad_9,tol_Res,tol_stop)	
tid = time()-tid;
println("OD ", tid);
Uh = ϕ'*UGS;
save("./Results/Rough/J_method_N_"*string(N)*"_Nh"*string(Nh)*"_l"*string(ℓ)*".jld","UGS",UGS,"Ẽ",Ẽ,"conv",conv,"tid",tid, "E", E,"Uh", Uh,"tid_assembly",tid_assembly)

end
end
