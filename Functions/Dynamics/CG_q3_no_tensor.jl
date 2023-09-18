function CG_q3(U0,T,Nt,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,TOL,max_it,β,ϵ,save_it,save_here)
# Lagrage FEM of order p for spatial discretization
k = T/Nt;
## numerical parameter


# cG(q)-scheme parameter
q = 3; #degree in time


τ,w = get_gauss_nodes(3,0,1); #Gauss integration
τ_t,w_t = get_gauss_nodes(2*3,0,1); # Gauss integration for nonlinear term

# Lagrange polynomials
li = [ x-> ((x-τ[2]).*(x-τ[3]))./((τ[1]-τ[2]).*(τ[1]-τ[3])); 
       x-> ((x-τ[1]).*(x-τ[3]))./((τ[2]-τ[1]).*(τ[2]-τ[3])); 
       x-> ((x-τ[1]).*(x-τ[2]))./((τ[3]-τ[1]).*(τ[3]-τ[2])) ];
   
li_hat = [ x-> -((x-τ[1]).*(x-τ[2]).*(x-τ[3]))./(τ[1].*τ[2].*τ[3]); 
           x-> (x.*(x-τ[2]).*(x-τ[3]))./(τ[1].*(τ[1]-τ[2]).*(τ[1]-τ[3])); 
           x-> (x.*(x-τ[1]).*(x-τ[3]))./(τ[2].*(τ[2]-τ[1]).*(τ[2]-τ[3]));  
           x-> (x.*(x-τ[1]).*(x-τ[2]))./(τ[3].*(τ[3]-τ[1]).*(τ[3]-τ[2])) ];
       
dli = [ x-> -(τ[2]+τ[3]-2*x)./((τ[1]-τ[2])*(τ[1]-τ[3])); 
        x-> -(τ[1]+τ[3]-2*x)./((τ[2]-τ[1])*(τ[2]-τ[3])); 
        x-> -(τ[1]+τ[2]-2*x)./((τ[3]-τ[1])*(τ[3]-τ[2])) ];

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
LU_Matrix_1 = lu(φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ ) ;
LU_Matrix_2 = lu(φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ );
LU_Matrix_3 = lu(φᵢφⱼ + 1im*ϵ*k*γ[3]*∇φᵢ∇φⱼ + 1im*k*γ[3]*Vφᵢφⱼ );

# time integration
# set intial starting point for fixed point iteration
ξⁱ_1 = zeros(SpaceDim_C)*0im;
ξⁱ_2 = zeros(SpaceDim_C)*0im;
ξⁱ_3 = zeros(SpaceDim_C)*0im;

U = zeros(SpaceDim_C)*0im;
copyto!(U,U0)

φᵢφⱼ_lu = lu(φᵢφⱼ);
tid = time();




for n = 1:Nt
    
    # fixed-point itreration
    Δ = 1;
    it = 0;

    copyto!(ξⁱ_1,U);
    copyto!(ξⁱ_2,U);
    copyto!(ξⁱ_3,U);
 tid = time();   
    while Δ > TOL && it<max_it
  
#  tid1 = time();      
        Û = c0[1]*U+c[1,1]*ξⁱ_1+c[2,1]*ξⁱ_2+c[3,1]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ;
        f₁ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
 
        
        Û = c0[2]*U+c[1,2]*ξⁱ_1+c[2,2]*ξⁱ_2+c[3,2]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₂ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
        
        Û = c0[3]*U+c[1,3]*ξⁱ_1+c[2,3]*ξⁱ_2+c[3,3]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₃ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
        
        Û = c0[4]*U+c[1,4]*ξⁱ_1+c[2,4]*ξⁱ_2+c[3,4]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
	ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₄ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
  
      
        Û = c0[5]*U+c[1,5]*ξⁱ_1+c[2,5]*ξⁱ_2+c[3,5]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₅ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
  
        Û = c0[6]*U+c[1,6]*ξⁱ_1+c[2,6]*ξⁱ_2+c[3,6]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₆ = ϕ*Assemble.ρuvᵢ(mesh,Quad8,Ûₕ,ρₕ);
 # println("TIME ASSEMBLY ", time()-tid1)
  #tid2 = time();
        # compute right hand side
        F1	= a[1]*φᵢφⱼ*U - 1im*k*β*( b[1,1]*f₁ + b[1,2]*f₂ + b[1,3]*f₃ +b[1,4]*f₄+b[1,5]*f₅+b[1,6]*f₆ );                                       
        F2	= a[2]*φᵢφⱼ*U - 1im*k*β*( b[2,1]*f₁ + b[2,2]*f₂ + b[2,3]*f₃ +b[2,4]*f₄+b[2,5]*f₅+b[2,6]*f₆ );
        F3 	= a[3]*φᵢφⱼ*U - 1im*k*β*( b[3,1]*f₁ + b[3,2]*f₂ + b[3,3]*f₃ +b[3,4]*f₄+b[3,5]*f₅+b[3,6]*f₆ );
 	
 	ξⁱ⁺¹_1 = LU_Matrix_1\F1;
        ξⁱ⁺¹_2 = LU_Matrix_2\F2;
        ξⁱ⁺¹_3 = LU_Matrix_3\F3;
        # update iteration paramter
 #  println("TIME SOLVE ", time()-tid2);     
          
        ΔU_1 = ξⁱ_1 - ξⁱ⁺¹_1; ΔU_2 = ξⁱ_2 - ξⁱ⁺¹_2; ΔU_3 = ξⁱ_3 - ξⁱ⁺¹_3
        Δ = sqrt(max( real(dot(ΔU_1,φᵢφⱼ*ΔU_1)),real(dot(ΔU_2,φᵢφⱼ*ΔU_2)),real(dot(ΔU_3,φᵢφⱼ*ΔU_3))) ) 
   
        copyto!(ξⁱ_1,ξⁱ⁺¹_1);
        copyto!(ξⁱ_2,ξⁱ⁺¹_2);
        copyto!(ξⁱ_3,ξⁱ⁺¹_3);
   
        
        it+=1;

    end
   println(time()-tid, " it ", it)
  
    # set result and calculate U(n+1)
   # 
    ξ_transform_1 = S⁻¹[1,1]*ξⁱ_1+S⁻¹[1,2]*ξⁱ_2+S⁻¹[1,3]*ξⁱ_3;
    ξ_transform_2 = S⁻¹[2,1]*ξⁱ_1+S⁻¹[2,2]*ξⁱ_2+S⁻¹[2,3]*ξⁱ_3;
    ξ_transform_3 = S⁻¹[3,1]*ξⁱ_1+S⁻¹[3,2]*ξⁱ_2+S⁻¹[3,3]*ξⁱ_3;
    
    Ua = li_hat[1](1/4)*U + li_hat[2](1/4)*ξ_transform_1+ li_hat[3](1/4)*ξ_transform_2+li_hat[4](1/4)*ξ_transform_3;   
    Ub = li_hat[1](1/2)*U + li_hat[2](1/2)*ξ_transform_1+ li_hat[3](1/2)*ξ_transform_2+li_hat[4](1/2)*ξ_transform_3; 
    Uc = li_hat[1](3/4)*U + li_hat[2](3/4)*ξ_transform_1+ li_hat[3](3/4)*ξ_transform_2+li_hat[4](3/4)*ξ_transform_3; 

    U = li_hat[1](1)*U + li_hat[2](1)*ξ_transform_1+ li_hat[3](1)*ξ_transform_2+li_hat[4](1)*ξ_transform_3;
    
          if(mod(n,save_it)==0)
		
		Uₕ = ϕ'*U;
		ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad8,Uₕ))
			 	
	 	E_nl = β/2*dot(ρ,φᵢφⱼ*ρ);
	 	#println(E_nl);
		Energy = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+E_nl;

		println("ENERGY"," ", (round((Energy)*10^7)/10^7).re,  " time ", n*k);
		
		save(save_here*"sol_SLOD_"*string(n)*".jld","U",U,"Ua",Ua,"Ub",Ub,"Uc",Uc);
		 
	end

   
end
return U
end

