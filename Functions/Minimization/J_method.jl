include("./GS.jl")


#----------------------------------------------------------------------------#





function J_METHOD_ωₕ(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,U,ϕ,max_it,φᵢφⱼ_chol,ωₕ,ω̃ₕ,ϵ,β,Quad,tol_shift,tol_stop)


ϕᵀ = sparse(ϕ');

E_exact = 0.;

conv_history = zeros(max_it,3);
Aᵥ = (ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ);

SpaceDim_C = size(ϕ,1); #coarse
SpaceDim_f = size(ϕ,2); #fine 	

Mass = sqrt(dot(U,(φᵢφⱼ*U)));
U = U/Mass
	
Uₕ = ϕᵀ*U;
ρ =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));

dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
U_interim = zeros(size(mesh.p,2))
U_interim[dofs_f] .= ϕ'*U;

Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
    	
Energy = U'*Aᵥ*U+β/2*ρ'*(φᵢφⱼ*ρ);
println("ENERGY START"," ", round((Energy)*10^7)/10^7, " Eᵧ ", Eᵧ, " Eb ",ρ'*(φᵢφⱼ*ρ));

tid_online = time();
	
TRI = sparse(LowerTriangular(vᵢvⱼ));

G1 = zeros(SpaceDim_C);
G2 = zeros(SpaceDim_C);
for n = 1:max_it

	tid = time();
	RHS = φᵢφⱼ*U;

        Uₕ = ϕ'*U;
        ρ = φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
        ρₕ = ϕ'*ρ;
	ρvᵢvⱼ = Assemble.ρφᵢφⱼ(ρₕ,ωₕ,ω̃ₕ,TRI);
        ρφᵢφⱼ = ϕ*ρvᵢvⱼ*ϕᵀ;

	# J = Aᵥ + 3*β*ρφᵢφⱼ-xy^T with 
	# x = ρφᵢφⱼ*U  and  y = 2*β*U'*φᵢφⱼ;
	J_sparse = Aᵥ + 3*β*ρφᵢφⱼ;  
	x = ρφᵢφⱼ*U;
	y = 2*β*U'*φᵢφⱼ;
		
	# Approximation of the eigenvalue and residual
	λ = U'*Aᵥ*U+β*(U'*ρφᵢφⱼ*U);
	Res = J_sparse*U - (y*U)*x - λ*(φᵢφⱼ*U);
	#Res = sqrt(dot(Res,∇φᵢ∇φⱼ\Res));
	Res = sqrt(dot(Res,Res));
		
	#println("Residual ", Res)
	# Shifting
	if Res < tol_shift
		# Rayleigh shifted
		σ = -λ;
	else
		σ = 0;
	end
		
		
	# Sherman–Morrison formula: (A-xy^T)^-1 = A^-1 +(A^-1xy^TA^-1)/(1-y^TA^-1x)
	# A = J_sparse + σ*φᵢφⱼ
	# x = ρφᵢφⱼ*U;
	# y = 2*β*U'*φᵢφⱼ;
	
	Jσ = J_sparse + σ*φᵢφⱼ;
	G1 = cg(Jσ,RHS);
	G2 = cg(Jσ,x);
	
	G = G1 + ((y*G1)/(1-y*G2))*G2;
	G = 1/(U'*(φᵢφⱼ*G))*G;
		
	if Res < tol_shift
		τ = 1;
	else

              Gₕ = ϕ'*G;
                ρ_G =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Gₕ,ωₕ));
                ρ_Gₕ = ϕ'*ρ_G;
		Gρ_Gvᵢ =ϕ*Assemble.u²uvᵢ(ωₕ,Gₕ,ρ_Gₕ);
		

		a0 = dot(U,Aᵥ*U);
		a1 = 2.0*dot(U,Aᵥ*G);
		a2 = dot(G,Aᵥ*G);

  
		b0 = β/2*dot(U,ρφᵢφⱼ*U);
		b1 = 2*β*dot(G,ρφᵢφⱼ*U);
		b2 = 3*β*dot(G,ρφᵢφⱼ*G);
		b3 = 2*β*dot(U,Gρ_Gvᵢ);
		b4 = β/2*dot(G,Gρ_Gvᵢ);
			
		z0 =   dot(U,φᵢφⱼ*U);
		z1 = 2*dot(U,φᵢφⱼ*G);
		z2 =   dot(G,φᵢφⱼ*G);
		
		# GOLDEN SECTION SEARCH
		
		τ,opt = Golden_section(a0,a1,a2,b0,b1,b2,b3,b4,z0,z1,z2);
	end
	

		
	U = (1-τ)*U+τ*G;
	U ./= sqrt(dot(U,(φᵢφⱼ*U)));
	
        Uₕ = ϕ'*U;
        ρ =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));

	out = U'*Aᵥ*U+β/2*(ρ'*(φᵢφⱼ*ρ));
		
	diff_E = Energy - out;	
	Energy = out;
	println("Energy: ", Energy, ", Residual: ", Res, ", λ: ", λ, " time ", round(100*(time()-tid))/100);
	conv_history[n,1:3] = [real(Energy), Res, λ]; 
	if(abs(diff_E)<10^-10);println( "time online ", time()-tid_online); 
		println("Residual ", Res)
			
		dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
		U_interim = zeros(size(mesh.p,2))
		U_interim[dofs_f] .= ϕ'*U;

		Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
		E_exact = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+β/2*Eᵧ;
		println(E_exact, " <---- Exact energy")
		break; 
	end
		
end
	
return U,Energy,E_exact,conv_history
end

function J_METHOD_Rot_ωₕ(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,φᵢLzφⱼ,U,ϕ,max_it,φᵢφⱼ_chol,ωₕ,ω̃ₕ,ϵ,Ω,β,Quad,tol_Switch,tol_Res,E_TOL)


ϕᵀ = sparse(ϕ');

conv_history = zeros(max_it,3);
φᵢLzφⱼ = real(φᵢLzφⱼ); #
Aᵥ = vcat(hcat(ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ, -Ω*φᵢLzφⱼ),hcat(Ω*φᵢLzφⱼ,ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ));
L = ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ+1im*Ω*φᵢLzφⱼ;

SpaceDim_C = size(ϕ,1);
SpaceDim_f = size(ϕ,2);

Mass = sqrt(dot(U,(φᵢφⱼ*U)));
U = U/Mass;
U_vec = [real(U);imag(U)];

Uₕ = ϕᵀ*U;

dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);	

	
tid_online = time();
	 
TRI = sparse(LowerTriangular(sparse(vᵢvⱼ)));

φᵢφⱼ_Bd = blockdiag(φᵢφⱼ,φᵢφⱼ)

max_it_IT = 100*size(φᵢφⱼ,1)

G = zeros(SpaceDim_C)*0im;
G1 = zeros(2*SpaceDim_C);
G2 = zeros(2*SpaceDim_C);


Res = 1.
Energy = 10^6;
E_exact = 1.0
N_it = 0



for n = 1:max_it
tid = time();

	if( (Res < tol_Switch)   );#ONLY SPLIT INTO REAL AND IM IF NECESSARY
		U_vec = [real(U);imag(U)];
	
		RHS = [φᵢφⱼ*real(U);φᵢφⱼ*imag(U)];
		
		Uₕ = ϕᵀ*U;
		ρ_Re²  = φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(real(Uₕ),ωₕ));
		ρ_Im²  = φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(imag(Uₕ),ωₕ));
		ρ_ReIm = φᵢφⱼ_chol\(ϕ*Assemble.u²uvᵢ(ωₕ,real(Uₕ),imag(Uₕ)));
		ρ  = ρ_Re² + ρ_Im²;



		ρφᵢφⱼ_Re²= ϕ*Assemble.ρφᵢφⱼ(ϕᵀ*ρ_Re²,ωₕ,ω̃ₕ,TRI)*ϕᵀ; 
		ρφᵢφⱼ_Im² = ϕ*Assemble.ρφᵢφⱼ(ϕᵀ*ρ_Im²,ωₕ,ω̃ₕ,TRI)*ϕᵀ; 
		ρφᵢφⱼ_ReIm = ϕ*Assemble.ρφᵢφⱼ(ϕᵀ*ρ_ReIm,ωₕ,ω̃ₕ,TRI)*ϕᵀ; 
		
#=
		ρφᵢφⱼ_Re² = LinearMap(x_h-> ϕ*Assemble.u²uvᵢ(ωₕ,ρ_Re²_h,x_h),SpaceDim_C,SpaceDim_f); 
		ρφᵢφⱼ_Im² = LinearMap(x_h-> ϕ*Assemble.u²uvᵢ(ωₕ,ρ_Im²_h,x_h),SpaceDim_C,SpaceDim_f); =#		
			Mρ = ρφᵢφⱼ_Re²+ρφᵢφⱼ_Im²

	   	
		#J_sparse = Aᵥ + β*vcat(hcat(3*ρφᵢφⱼ_Re²+ρφᵢφⱼ_Im²,2*LIm1),			hcat(2*LIm2,ρφᵢφⱼ_Re²+3*ρφᵢφⱼ_Im²));  
		
		J_sparse = Aᵥ + β*vcat(hcat(3*ρφᵢφⱼ_Re²+ρφᵢφⱼ_Im²,2*ρφᵢφⱼ_ReIm),hcat(2*ρφᵢφⱼ_ReIm,ρφᵢφⱼ_Re²+3*ρφᵢφⱼ_Im²));  
		x = [Mρ*real(U);Mρ*imag(U)];
		y = 2*β*[φᵢφⱼ*real(U);φᵢφⱼ*imag(U)]';
	
		
		# Approximation of the eigenvalue and residual
		λ = dot(U_vec,J_sparse*U_vec)-dot(U_vec,x)*dot(y,U_vec);
		λ_Sobolev = dot(U_vec,Aᵥ*U_vec)+β*dot(U_vec,x);
		
		
		Res = J_sparse*U_vec - dot(y,U_vec)*x - λ*(φᵢφⱼ_Bd*U_vec);
		Res = sqrt(dot(Res,Res));
	
		# Shifting:
		#σ = -λ
		Jσ = (J_sparse-λ*φᵢφⱼ_Bd)
		
		pl = DiagonalPreconditioner(Jσ)
		idrs!(G1,Jσ,RHS,Pl = pl,reltol = 1e-8, maxiter = SpaceDim_C*10)
		idrs!(G2,Jσ,x,Pl = pl, reltol = 1e-8, maxiter = SpaceDim_C*10)
	
		G1+= dot(G1,y)/(1-dot(y,G2))*G2;   
		G .= G1[1:SpaceDim_C] + 1im*G1[(1+SpaceDim_C):2*SpaceDim_C];

		tol_Switch/=1.1;
		N_it +=1;
		
		U = G;
		U ./= sqrt(real(dot(U,φᵢφⱼ*U)));
		
		ρ =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
 		
		opt = real(U'*(L*U)+β/2*ρ'*(φᵢφⱼ*ρ));
	else 
		N_it +=1;
		RHS = φᵢφⱼ*U;
	
		Uₕ = ϕᵀ*U;
		ρ =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
 		ρₕ =ϕᵀ*ρ; #consider changing
 	
		ρvᵢvⱼ = Assemble.ρφᵢφⱼ(ρₕ,ωₕ,ω̃ₕ,TRI);
        	Mρ = ϕ*(ρvᵢvⱼ*ϕᵀ);
		#Mρ =LinearMap(x-> ϕ*(ρvᵢvⱼ*(ϕᵀ*x)),SpaceDim_C,SpaceDim_C);
		λ = real(U'*(L*U)+β*U'*(Mρ*U));
		d = (L*U)+β*(Mρ*U)-λ*(φᵢφⱼ*U); Res = sqrt(real(dot(d,d)));
		
		Mit = L+β*Mρ;
		pl = DiagonalPreconditioner(Mit);
		idrs!(G,Mit,RHS,reltol = 10^-7,Pl=pl); 
		#idrs!(G,L+β*Mρ,RHS,reltol = 10^-7); 
		
		Gₕ = ϕᵀ*G;
		ρ_G =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Gₕ,ωₕ));
		ρ_Gₕ = ϕᵀ*ρ_G;
		Gρ_Gvᵢ =ϕ*Assemble.u²uvᵢ(ωₕ,Gₕ,ρ_Gₕ);
		

		a0 = real(dot(U,L*U));
		a1 = 2.0*real(dot(U,L*G));
		a2 = real(dot(G,L*G));

  
  		b0 = β/2*real(dot(U,Mρ*U));
		b1 = 2*β*real(dot(G,Mρ*U));
		GU = φᵢφⱼ_chol\(ϕ*real(Assemble.u²uvᵢ(ωₕ,ϕᵀ*U,conj(Gₕ))));
		b2 = β*real(dot(G,Mρ*G)+2*dot(GU,φᵢφⱼ*GU));

		b3 = 2*β*real(dot(U,Gρ_Gvᵢ));
		b4 = β/2*real(dot(G,Gρ_Gvᵢ));
		
		z0	= real(dot(U,φᵢφⱼ*U));
		z1 	= 2*real(dot(U,φᵢφⱼ*G));
		z2	= real(dot(G,φᵢφⱼ*G));
	
	
	
		#----------------GOLDEN SECTION SEARCH WITH ASSUMPTIONS ON PROXIMITY TO ROOT-----------------

		τ,opt = Golden_section(a0,a1,a2,b0,b1,b2,b3,b4,z0,z1,z2)
		
		U = (1-τ)*U+τ*G;
		U ./= sqrt(real(dot(U,φᵢφⱼ*U)));
	
	end
	
	

	diff_E = Energy - opt;	
	Energy = opt;
	println(n, ": Energy ", real(Energy), " λ ", λ, " Res ", Res , " time ", round(100*(time()-tid))/100)
	
	conv_history[n,:] = [real(Energy), λ, Res]; 
	if( (abs(Res)<tol_Res) | (N_it == max_it) | ( abs(diff_E)<E_TOL) );
		println( "time online ", time()-tid_online); 

		dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
		U_interim = 0im*zeros(size(mesh.p,2))
		U_interim[dofs_f] .= ϕ'*U;

		Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
		E_exact = real(U'*((ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ+1im*Ω*φᵢLzφⱼ)*U)+β/2*Eᵧ);
		
		N_it = n;
		break; 
	end
end

return U,Energy,E_exact,conv_history[1:N_it,:]
end






function J_METHOD_ωₕ_IterativeSolvers(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,U,ϕ,max_it,ωₕ,ω̃ₕ,ϵ,β,Quad,tol_Res,tol_stop)

ϕᵀ = sparse(ϕ');

E_exact = 0.;

conv_history = zeros(max_it,3);
Aᵥ = (ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ);


SpaceDim_C = size(ϕ,1); #coarse
SpaceDim_f = size(ϕ,2); #fine 	

Mass = sqrt(dot(U,(φᵢφⱼ*U)));
U = U/Mass
	
Uₕ = ϕ'*U;
ρ =cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));

dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
U_interim = zeros(size(mesh.p,2))
U_interim[dofs_f] .= ϕ'*U;

Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
    	
Energy = U'*Aᵥ*U+β/2*ρ'*(φᵢφⱼ*ρ);
println("ENERGY START"," ", round((Energy)*10^7)/10^7, " Eᵧ ", Eᵧ, " Eb ",ρ'*(φᵢφⱼ*ρ));

tid_online = time();
	
TRI = sparse(LowerTriangular(vᵢvⱼ));
pr = CholeskyPreconditioner(φᵢφⱼ, 2)

for n = 1:max_it

	tid = time();
	RHS = φᵢφⱼ*U;

        Uₕ = ϕ'*U;
       # ρ = φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
        ρ = cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Uₕ,ωₕ),Pl=pr);
        
        ρₕ = ϕ'*ρ;
	ρvᵢvⱼ = Assemble.ρφᵢφⱼ(ρₕ,ωₕ,ω̃ₕ,TRI);
        ρφᵢφⱼ = ϕ*ρvᵢvⱼ*ϕᵀ;

	# J = Aᵥ + 3*β*ρφᵢφⱼ-xy^T with 
	# x = ρφᵢφⱼ*U  and  y = 2*β*U'*φᵢφⱼ;
	J_sparse = Aᵥ + 3*β*ρφᵢφⱼ;  
	x = ρφᵢφⱼ*U;
	y = 2*β*U'*φᵢφⱼ;
		
	# Approximation of the eigenvalue and residual
	λ = U'*Aᵥ*U+β*(U'*ρφᵢφⱼ*U);
	Res = J_sparse*U - (y*U)*x - λ*(φᵢφⱼ*U);
	#Res = sqrt(dot(Res,∇φᵢ∇φⱼ\Res));
	Res = sqrt(dot(Res,Res));
		
	#println("Residual ", Res)
	# Shifting
	if Res < tol_Res
		# Rayleigh shifted
		σ = -λ;
	else
		σ = 0;
	end
		
		
	# Sherman–Morrison formula: (A-xy^T)^-1 = A^-1 +(A^-1xy^TA^-1)/(1-y^TA^-1x)
	# A = J_sparse + σ*φᵢφⱼ
	# x = ρφᵢφⱼ*U;
	# y = 2*β*U'*φᵢφⱼ;
	Jσ = J_sparse + σ*φᵢφⱼ;	
	p = DiagonalPreconditioner(Jσ);

	G1 = cg(Jσ,RHS,Pl=p); #
	G2 = cg(Jσ,x,Pl=p);
	G = G1 + ((y*G1)/(1-y*G2))*G2;
	G = 1/(U'*(φᵢφⱼ*G))*G;
		
	if Res < tol_Res
		τ = 1;
	else

              Gₕ = ϕ'*G;
               # ρ_G =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Gₕ,ωₕ));
                ρ_G =cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Gₕ,ωₕ),Pl=pr);
                ρ_Gₕ = ϕ'*ρ_G;
		Gρ_Gvᵢ =ϕ*Assemble.u²uvᵢ(ωₕ,Gₕ,ρ_Gₕ);
		

		a0 = dot(U,Aᵥ*U);
		a1 = 2.0*dot(U,Aᵥ*G);
		a2 = dot(G,Aᵥ*G);

  
		b0 = β/2*dot(U,ρφᵢφⱼ*U);
		b1 = 2*β*dot(G,ρφᵢφⱼ*U);
		b2 = 3*β*dot(G,ρφᵢφⱼ*G);
		b3 = 2*β*dot(U,Gρ_Gvᵢ);
		b4 = β/2*dot(G,Gρ_Gvᵢ);
			
		z0 =   dot(U,φᵢφⱼ*U);
		z1 = 2*dot(U,φᵢφⱼ*G);
		z2 =   dot(G,φᵢφⱼ*G);
		
		# GOLDEN SECTION SEARCH
		
		τ,opt = Golden_section(a0,a1,a2,b0,b1,b2,b3,b4,z0,z1,z2);
	end
		
	U = (1-τ)*U+τ*G;
	U ./= sqrt(dot(U,(φᵢφⱼ*U)));
	
        Uₕ = ϕ'*U;
       # ρ =φᵢφⱼ_chol\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
	 ρ =cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Uₕ,ωₕ),Pl=pr);

	out = U'*Aᵥ*U+β/2*(ρ'*(φᵢφⱼ*ρ));
		
	diff_E = Energy - out;	
	Energy = out;
	println("Energy: ", Energy, ", Residual: ", Res, ", λ: ", λ);
	conv_history[n,1:3] = [real(Energy), Res, λ]; 
	
	if(abs(Res)<tol_stop);println( "time online ", time()-tid_online); 
		println("Residual ", Res)
			
		dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
		U_interim = zeros(size(mesh.p,2))
		U_interim[dofs_f] .= ϕ'*U;

		Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
		E_exact = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+β/2*Eᵧ;
		println(E_exact, " <---- Exact energy")
		break; 
	end
	#println(time()-tid);	
end
	
return U,Energy,E_exact,conv_history
end






