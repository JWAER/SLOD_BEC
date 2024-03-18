function CG_q2_ωₕ_IterativeSolvers(U0,T,Nt,ϕ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,ω,TOL,max_it,β,ϵ,save_it,save_here,Quad)
# Lagrage FEM of order p for spatial discretization
k = T/Nt;
## numerical parameter


# cG(q)-scheme parameter
q = 2; #degree in time


τ,w = get_gauss_nodes(q,0,1); #Gauss integration
τ_t,w_t = get_gauss_nodes(2*q,0,1); # Gauss integration for nonlinear term

# Lagrange polynomials
li,li_hat,dli = lagrange_basis(q,τ);

# compute RK coefficients
m = zeros(q,q);
for i = 1:q
    for j = 1:q
        if i==j
            m[i,j] = (w[i]/τ[j])*(1 + τ[i]*dli[j](τ[i]));
        else
            m[i,j] = (w[i]/τ[j])*τ[i]*dli[j](τ[i]);
        end
    end
end


m⁻¹ = m\I;
γ,S⁻¹ = eigen(m⁻¹*diagm(w));
D = diagm(γ);
S = S⁻¹\I;

# a coefficients
a = zeros(q,1)*0im;
for i = 1:q
    a[i] = sum(S[i,:]);
end

# b coeeficitents
b = zeros(q,2*q)*0im;
SMinv = S*m⁻¹;
for i = 1:q
    for m = 1:2*q
        for j = 1:q
            b[i,m] = b[i,m] + SMinv[i,j]*li[j](τ_t[m]);
        end
        b[i,m] = w_t[m]*b[i,m];
    end
end
# c coefficients
c = zeros(q,2*q)*0im;
for j = 1:q
    for m = 1:2*q
        for i = 1:q
            c[j,m] = c[j,m] + li_hat[i+1](τ_t[m])*S⁻¹[i,j];
        end
    end
end
c0 = zeros(2*q,1);
for m = 1:2*q
   c0[m] = li_hat[1](τ_t[m]); 
end


# compute LHS and LU-decomposition
M1=φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ ;
M2=φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ ;

M1_precond = ilu(M1,τ = 0.00001);
M2_precond = ilu(M2,τ = 0.00001);


# time integration
# set intial starting point for fixed point iteration
ξⁱ_1 = zeros(SpaceDim_C)*0im;
ξⁱ_2 = zeros(SpaceDim_C)*0im;

U = zeros(SpaceDim_C)*0im;
copyto!(U,U0)

tid = time();


pr = CholeskyPreconditioner(φᵢφⱼ, 2)

for n = 1:Nt
    
    # fixed-point itreration
    Δ = 1;
    it = 0;

    copyto!(ξⁱ_1,U);
    copyto!(ξⁱ_2,U);


tiiid = time();
    while Δ > TOL && it<max_it
  
        Û = c0[1]*U+c[1,1]*ξⁱ_1+c[2,1]*ξⁱ_2;
        Ûₕ = ϕ'*Û;
        ρ = cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Ûₕ,ωₕ),Pl=pr) 
        ρₕ = ϕ'*ρ;
        f₁ = ϕ*Assemble.u²uvᵢ(ωₕ,Ûₕ,ρₕ)
 
        
        Û = c0[2]*U+c[1,2]*ξⁱ_1+c[2,2]*ξⁱ_2;
        Ûₕ = ϕ'*Û;
        ρ = cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Ûₕ,ωₕ),Pl=pr) 
        ρₕ = ϕ'*ρ; 
        f₂ = ϕ*Assemble.u²uvᵢ(ωₕ,Ûₕ,ρₕ)
        
        Û = c0[3]*U+c[1,3]*ξⁱ_1+c[2,3]*ξⁱ_2;
        Ûₕ = ϕ'*Û;
        ρ =  cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Ûₕ,ωₕ),Pl=pr)
        ρₕ = ϕ'*ρ; 
        f₃ = ϕ*Assemble.u²uvᵢ(ωₕ,Ûₕ,ρₕ)
        
        Û = c0[4]*U+c[1,4]*ξⁱ_1+c[2,4]*ξⁱ_2;
        Ûₕ = ϕ'*Û;
	ρ =  cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Ûₕ,ωₕ),Pl=pr)
        ρₕ = ϕ'*ρ; 
        f₄ = ϕ*Assemble.u²uvᵢ(ωₕ,Ûₕ,ρₕ)
        # compute right hand side
        F1	= a[1]*φᵢφⱼ*U - 1im*k*β*( b[1,1]*f₁ + b[1,2]*f₂ + b[1,3]*f₃ +b[1,4]*f₄);                                       
        F2	= a[2]*φᵢφⱼ*U - 1im*k*β*( b[2,1]*f₁ + b[2,2]*f₂ + b[2,3]*f₃ +b[2,4]*f₄);
     
        ξⁱ⁺¹_1 = idrs(M1,F1,Pl=M1_precond,reltol=10^-11);
        ξⁱ⁺¹_2 = idrs(M2,F2,Pl=M2_precond,reltol=10^-11);
        
        
        # update iteration paramter
          
        ΔU_1 = ξⁱ_1 - ξⁱ⁺¹_1; ΔU_2 = ξⁱ_2 - ξⁱ⁺¹_2; 
        Δ = sqrt(max( real(dot(ΔU_1,φᵢφⱼ*ΔU_1)),real(dot(ΔU_2,φᵢφⱼ*ΔU_2))) ) 
   
        copyto!(ξⁱ_1,ξⁱ⁺¹_1);
        copyto!(ξⁱ_2,ξⁱ⁺¹_2);
   
        
        it+=1;

    end
  println("one time step ", time()-tiiid, " it ", it)
  
    # set result and calculate U(n+1)
   # 
   tid = time();
    ξ_transform_1 = S⁻¹[1,1]*ξⁱ_1+S⁻¹[1,2]*ξⁱ_2;
    ξ_transform_2 = S⁻¹[2,1]*ξⁱ_1+S⁻¹[2,2]*ξⁱ_2;
      
    

    U = li_hat[1](1)*U + li_hat[2](1)*ξ_transform_1+ li_hat[3](1)*ξ_transform_2;
    
          if(mod(n,save_it)==0)
		
		Uₕ = ϕ'*U;
		ρ = cg(φᵢφⱼ,ϕ*Assemble.ρvᵢ(Uₕ,ωₕ),Pl=pr)
			 	
	 	E_nl = β/2*dot(ρ,φᵢφⱼ*ρ);
	 	#println(E_nl);
		Energy = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+E_nl;

		U_interim = 0im*zeros(size(mesh.p,2))
		U_interim[mesh.dofs] .= Uₕ;

		#Eᵧ =  Assemble.NL_Energy(U_interim,Quad,mesh);
		#E_exact = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+β/2*Eᵧ;

		println("ENERGY"," ", (round((Energy)*10^7)/10^7).re,  " time ", n*k);
#		println("ENERGY ", E_exact)
		
		save(save_here*"sol_SLOD_"*string(n)*".jld","U",U,"Energy",Energy)
		 
	end
println("all else ", time()-tid);


   
end
return U
end



