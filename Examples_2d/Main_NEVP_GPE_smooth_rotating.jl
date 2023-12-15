include("../Functions/Dependencies_2d.jl");

#-----------------------------define constants and potentials-----------------------------------#

β = 1000
Ω = 0.85
N = 80
ℓ = 2;
Nh = 1; #refinement, for smooth problems Nh = 1;
box = 20; #will give a square domain from -box/2 to box/2
ϵ = 1.0;
γ_x = 1.0;
γ_y = 1.0;
Vs(x) = 1/2*dot(x,x);

dynamics = !true
Vₜ(x) = 1/2*dot(x,x);#potential in time-dependent problem
	

#--------------------------Generate mesh ------------------------------------------
	mesh = SLOD_MESH.Generate(N,Nh,box,ℓ);
	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);

#------------------------------Assemble Matrix Quantities------------------------------------#
tid_assembly = time();
	∇vᵢ∇vⱼ ,vᵢvⱼ,vᵢvⱼV,P,vᵢLzvⱼ = Assemble.All(mesh,Quad_4,Quad_6,Quad_9,Vs,Ω)

	tid = time();
	ϕ = (Assemble.ϕ_canonical(mesh,∇vᵢ∇vⱼ,vᵢvⱼ,P,ℓ)[:,dofs_f]);
	println("Computed SLOD space in ", time()-tid);
	SpaceDim_C = size(ϕ,1);
	
	
#---------------------------------------- Compute Tensor ----------------------------------------#

	
	tid = time();
	ωₕ=ω_module.PreAllocateW(vᵢvⱼ[dofs_f,dofs_f])
	ω_module.Compute_Wh(mesh,ωₕ);

	ω̃ₕ  = ω_module.Compute_tilde(ωₕ)
	

#------------------------- Restrict to H^1_0 and compute SLOD-versions-------------------------------

	∇vᵢ∇vⱼ = ∇vᵢ∇vⱼ[dofs_f,dofs_f]
	vᵢvⱼ = vᵢvⱼ[dofs_f,dofs_f];
	vᵢvⱼV =vᵢvⱼV[dofs_f,dofs_f]
	vᵢLzvⱼ =vᵢLzvⱼ[dofs_f,dofs_f];
	

	∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ*ϕ')
	φᵢφⱼ = ϕ*(vᵢvⱼ*ϕ')
	Vφᵢφⱼ = ϕ*(vᵢvⱼV*ϕ');
	φᵢLzφⱼ = ϕ*(vᵢLzvⱼ*ϕ');
	

	tid_assembly = time()-tid_assembly;

#-----------------------------------Solve for minimizier --------------------------------------#

#---------------------------------Initial guess -----------------------------------------------#

	φᵢφⱼ_lu = lu(φᵢφⱼ); # compute lu once


	UGS = Assemble.initial_guess(Ω,mesh,φᵢφⱼ_lu,vᵢvⱼ,dofs_f,ϕ)


#----------------------------------------------------------------------------------------------#
	tid = time()
	max_it = 10000
	tol_stop = 1e-6; 
	tol_shift = 3e-3;
	tol_E = 1e-8;
	UGS, E,conv = J_METHOD_Rot_ωₕ(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,φᵢLzφⱼ,UGS,ϕ,max_it,φᵢφⱼ_lu,ωₕ,ω̃ₕ ,ϵ,Ω,β,Quad_9,tol_shift,tol_stop,tol_E)
	tid = time()-tid;
	Uₕ = ϕ'*UGS;
#	
	save("./Results/Smooth/J_method_rotating_N_"*string(N)*"_Omega_"*string(Ω)*"_beta_"*string(β)*".jld","UGS",UGS,"conv",conv,"tid",tid, "E", E,"Uₕ", Uₕ,"tid_assembly",tid_assembly)
#--------------------- Solve time-dependent problem ----------------------------------

ω̃ₕ  = 0;
	if(dynamics)
	
		Vₜvᵢvⱼ = Assemble.Vvᵢvⱼ(Vₜ,mesh,Quad_9);
		Vₜφᵢφⱼ = ϕ*(Vₜvᵢvⱼ[dofs_f,dofs_f]*ϕ');
		println("Computed smooth part in ", time()-tid);


		T = 1; Nt = 2^8; max_it = 20; save_it = 2; save_here = "./Results/"; TOL = 10^-8;
		CG_q3_ωₕ(UGS,T,Nt,∇φᵢ∇φⱼ,φᵢφⱼ,Vₜφᵢφⱼ,ωₕ,TOL,max_it,β,ϵ,save_it,save_here)
	end


#-----------------------Plotting using PlotlyJS --------------------------------------
using PlotlyJS


ρ = zeros(size(mesh.p,2));
for i = 1:length(dofs_f); ρ[dofs_f[i]]= abs(Uₕ[i])^2; end

ρ = reshape(ρ,3N+1,3N+1);

plot(heatmap(z=ρ,colorscale= "Jet"))





