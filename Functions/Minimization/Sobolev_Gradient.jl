include("./GS.jl")
function SOBOLEV_GRADIENT(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,U,ϕ,max_it,φᵢφⱼ_lu,ωₕ,ω̃ₕ,ϵ,β)

Conv_history = zeros(max_it);

Aᵥ = (ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ);

SpaceDim_C = size(ϕ,1);
SpaceDim_f = size(ϕ,2);	

	Uₕ = ϕ'*U;
	ρ =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
 
    	    	
    	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	U_interim = 0im*zeros(size(mesh.p,2))
	U_interim[dofs_f] .= ϕ'*U;

	ENL =  Assemble.NL_Energy(U_interim,Quad_9,mesh);
    	
    	

	Energy = U'*Aᵥ*U+β/2*ρ'*(φᵢφⱼ*ρ)
	println("ENERGY START"," ", round((Energy)*10^7)/10^7, " ENL ", ENL, " Eb ",ρ'*(φᵢφⱼ*ρ));



Mass = sqrt(dot(U,(φᵢφⱼ*U)));

U = U/Mass
#COMPUTE G(z) ---------------------
G = zeros(SpaceDim_C);


TRI = sparse(LowerTriangular(vᵢvⱼ));

for n = 1:max_it

	RHS = φᵢφⱼ*U;
	
	#Rho_U = @sync @distributed (+) for i = workers(); Assemble.ρvᵢ_parallel(U); end
	Uₕ = ϕ'*U;
	ρ =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
 	ρₕ =ϕ'*ρ; #consider changing
 	
	ρvᵢvⱼ = Assemble.ρφᵢφⱼ(ρₕ,ωₕ,ω̃ₕ,TRI);
        ρφᵢφⱼ = ϕ*(ρvᵢvⱼ*ϕ');

	ρuφᵢ = ρφᵢφⱼ*U; 
	cg!(G,Aᵥ+β*ρφᵢφⱼ,RHS,reltol = 10^-8); 
	Gₕ = ϕ'*G;
	
	Gρφᵢ = ϕ*(Assemble.u²uvᵢ(ωₕ,Gₕ,ρₕ)); #WEAK REPRESENTATION


	a0 = dot(U,Aᵥ*U);
	a1 = 2.0*dot(U,Aᵥ*G);
	a2 = dot(G,Aᵥ*G);

        ρ_G =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(Gₕ,ωₕ));
        ρ_Gₕ = ϕ'*ρ_G;
	Gρ_Gφᵢ =ϕ*Assemble.u²uvᵢ(ωₕ,Gₕ,ρ_Gₕ);
		
	b0 = β/2*dot(U,ρφᵢφⱼ*U);
	b1 = 2*β*dot(G,ρφᵢφⱼ*U);
	b2 = 3*β*dot(G,ρφᵢφⱼ*G);
	b3 = 2*β*dot(U,Gρ_Gφᵢ);
	b4 = β/2*dot(G,Gρ_Gφᵢ);
			
	z0 =   dot(U,φᵢφⱼ*U);
	z1 = 2*dot(U,φᵢφⱼ*G);
	z2 =   dot(G,φᵢφⱼ*G);

	τ,opt = Golden_section(a0,a1,a2,b0,b1,b2,b3,b4,z0,z1,z2)
	
	
	U = (1-τ )*U+τ *G
	U ./= sqrt(dot(U,φᵢφⱼ*U));



#------------------------------------------------------------------------------#

diff_E = Energy - opt	
Energy = opt;
Conv_history[n] = Energy;
println("Energy ", Energy);
if(abs(diff_E)<10^-10);#println( "time online ", time()-tid_online); 

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	U_interim = 0im*zeros(size(mesh.p,2))
	U_interim[dofs_f] .= ϕ'*U;

	ENL =  Assemble.NL_Energy(U_interim,Quad_9,mesh);
	E_exact = U'*(Aᵥ*U)+β/2*ENL;
	println(E_exact, " <---- Exact energy")
break; end

end

return U,Energy,Conv_history


end

















#----------------------------------------------------------------------------#


function SOBOLEV_GRADIENT_Rot(mesh,vᵢvⱼ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,MΩLOD,U,ϕ,max_it,φᵢφⱼ_lu,ωₕ,ω̃ₕ,ϵ,Ω,β)

conv_history = zeros(max_it,2);

Aᵥ = (ϵ*∇φᵢ∇φⱼ+Vφᵢφⱼ+1im*Ω*MΩLOD);

SpaceDim_C = size(ϕ,1);

	ρ =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(ϕ'U,ωₕ));
 
    	    	
    	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	U_interim = 0im*zeros(size(mesh.p,2))
	U_interim[dofs_f] .= ϕ'*U;

	ENL =  Assemble.NL_Energy(U_interim,Quad_9,mesh);

	Energy = U'*Aᵥ*U+β/2*ρ'*(φᵢφⱼ*ρ)
	println("ENERGY START"," ", round((Energy)*10^7)/10^7, " ENL ", ENL, " Eb ",ρ'*(φᵢφⱼ*ρ));



Mass = sqrt(dot(U,(φᵢφⱼ*U)));

U = U/Mass
#COMPUTE G(z) ---------------------
G = zeros(SpaceDim_C)*0im;


tid_online = time();


TRI = sparse(LowerTriangular(vᵢvⱼ));

for n = 1:max_it

	RHS = φᵢφⱼ*U;
	
	#Rho_U = @sync @distributed (+) for i = workers(); Assemble.ρvᵢ_parallel(U); end
	Uₕ = ϕ'*U;
	ρ =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(Uₕ,ωₕ));
 	ρₕ =ϕ'*ρ; #consider changing
 	
	ρvᵢvⱼ = Assemble.ρφᵢφⱼ(ρₕ,ωₕ,ω̃ₕ,TRI);
        ρφᵢφⱼ = ϕ*(ρvᵢvⱼ*ϕ');

    	λ = U'*Aᵥ*U+β*ρ'*(φᵢφⱼ*ρ);

	cg!(G,Aᵥ+β*ρφᵢφⱼ,RHS,reltol = 10^-8); #BADLY CONDITONED DUE TO BDRY..
	Gₕ = ϕ'*G;
	
	
	Gρφᵢ = ϕ*(Assemble.u²uvᵢ(ωₕ,Gₕ,ρₕ)); #WEAK REPRESENTATION

	
	a0 = real(dot(U,Aᵥ*U));
	a1 = real(2.0*dot(U,Aᵥ*G));
	a2 = real(dot(G,Aᵥ*G));

        ρ_G =φᵢφⱼ_lu\(ϕ*Assemble.ρvᵢ(Gₕ,ωₕ));
        ρ_Gₕ = ϕ'*ρ_G;
	Gρ_Gφᵢ =ϕ*Assemble.u²uvᵢ(ωₕ,Gₕ,ρ_Gₕ);
		
	b0 = β/2*real(dot(U,ρφᵢφⱼ*U));
	b1 = 2*β*real(dot(G,ρφᵢφⱼ*U));
	b2 = 3*β*real(dot(G,ρφᵢφⱼ*G));
	b3 = 2*β*real(dot(U,Gρ_Gφᵢ));
	b4 = β/2*real(dot(G,Gρ_Gφᵢ));
			
	z0 =   real(dot(U,φᵢφⱼ*U));
	z1 = 2*real(dot(U,φᵢφⱼ*G));
	z2 =   real(dot(G,φᵢφⱼ*G));

	τ,opt = Golden_section(a0,a1,a2,b0,b1,b2,b3,b4,z0,z1,z2)
	
	#println(τ)
	U = (1-τ )*U+τ *G
	U ./= sqrt(dot(U,φᵢφⱼ*U));



#------------------------------------------------------------------------------#

diff_E = Energy - opt	
Energy = opt;
println(n, ": Energy ", Energy, ", λ ", λ);
conv_history[n,:] = [real(Energy), real(λ)];
if(abs(diff_E)<10^-10);#println( "time online ", time()-tid_online); 

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	U_interim = 0im*zeros(size(mesh.p,2))
	U_interim[dofs_f] .= ϕ'*U;

	ENL =  Assemble.NL_Energy(U_interim,Quad_9,mesh);
	E_exact = U'*(Aᵥ*U)+β/2*ENL;
	println(E_exact, " <---- Exact energy")
break; end
#end
end

return U,Energy,conv_history
end

#-------------------------------------------------------------------------------------#

