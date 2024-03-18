include("../Functions/Dependencies_2d.jl");


#----------------------------Generate or load square domain -------------------------------------#
ℓ = 2;
Nh = 5; #refinement, for smooth problems Nh = 1;
Ns = [12 24 48 96 192];

for it = 1:5

βs = 100;
Ω = 0.0;
ϵ = 1/2; 
box = 12; #will give a square domain from -6 to 6
N = Ns[it]
Vd_refinement_factor = 6#Vd_refinement_factors[it];

mesh = SLOD_MESH.Generate(N,Nh,box,ℓ);

Vd(x) = 5 + min.(floor.(2*sin(π*x[1]/3)*sin(π*x[2]/3)),1);


γ_x = 1.0;
γ_y = 1.0;
Vs(x) = 1/2*dot(x,x);

dynamics = false;
βt = 20;
Vₜ(x) = 5 + min.(floor.(2*sin(π*x[1]/1)*sin(π*x[2]/1)),1);
#------------------------------Compute All Matrix Quantities------------------------------------#
#First compute in H¹
sparsity = Assemble.MatrixSparsity(mesh);

	tid = time();
	∇vᵢ∇vⱼ = Assemble.∇vᵢ∇vⱼ(mesh,Quad_4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad_6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	#vᵢvⱼVd = Assemble.Vvᵢvⱼ_piecewise_constant(Vd,mesh,Quad_9,Vd_refinement_factor,sparsity);
	vᵢvⱼVd = Assemble.Vvᵢvⱼ_example2(Vd,mesh,Quad_9,Vd_refinement_factor,sparsity);
	println("Computed rough matrix in ", time()-tid);

	tid = time();
	vᵢvⱼVs = Assemble.Vvᵢvⱼ(Vs,mesh,Quad_9,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.P(mesh)
	println("Computed P in ", time()-tid);



vᵢvⱼV = vᵢvⱼVs+vᵢvⱼVd;
Aᵥ = ϵ*∇vᵢ∇vⱼ+vᵢvⱼVd;

tid = time();
ϕ = Assemble.ϕ(mesh,Aᵥ,vᵢvⱼ,P,ℓ); #assembles EVERY SLOD-function ... possibly ok in 2d
println("Computed SLOD space in ", time()-tid);

#------------------------- Restrict to H^1_0   ----------------------------------
	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	
	ϕ = sparse(ϕ)[:,dofs_f];
	
	
	∇vᵢ∇vⱼ = ∇vᵢ∇vⱼ[dofs_f,dofs_f]
	vᵢvⱼ = vᵢvⱼ[dofs_f,dofs_f];
	vᵢvⱼV =vᵢvⱼV[dofs_f,dofs_f]
	

	∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ*ϕ')
	φᵢφⱼ = ϕ*(vᵢvⱼ*ϕ'); φᵢφⱼ +=φᵢφⱼ'; φᵢφⱼ /=2;
	Vφᵢφⱼ = ϕ*(vᵢvⱼV*ϕ');
	
	
#------------------------------------------------------------------------------------------------#

#---------------------------------------- Compute Tensor ----------------------------------------#

	tid = time();
	ωₕ=ω_module.PreAllocateW(vᵢvⱼ)
	ω_module.Compute_Wh(mesh,ωₕ);

	ω̃ₕ  = ω_module.Compute_tilde(ωₕ)


#----------------------------------------------------------------------------------------------#

#-----------------------------------Solve for minimizier --------------------------------------#


φᵢφⱼ_chol = cholesk(φᵢφⱼ); # compute lu once
UGS = Assemble.initial_guess(Ω,mesh,φᵢφⱼ_chol,vᵢvⱼ,dofs_f,ϕ)

#----------------------------------------------------------------------------------------------#
#max_it = 20
tid = time()
max_it = 1000
tol_shift = 0.05;
tol_stop = 1e-7;
UGS, E,conv,E_exact = J_METHOD_ωₕ(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,UGS,ϕ,max_it,φᵢφⱼ_chol,ωₕ,ω̃ₕ ,ϵ,βs,Quad_9,tol_shift,tol_stop)	  
tid = time()-tid;
println("OD ", tid);
Uh = ϕ'*UGS;
save("./Results/Rough/J_method_N_"*string(N)*"_Nh"*string(Nh)*"_l"*string(ℓ)*".jld","UGS",UGS,"E",E,"conv",conv,"tid",tid, "E_exact", E_exact,"Uh", Uh)

end

