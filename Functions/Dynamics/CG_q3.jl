function CG_q3(U0,T,Nt,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,ω,TOL,max_it,β,ϵ,save_it,save_here)
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


mi = m\I;
γ,Si = eigen(mi*diagm(w));
D = diagm(γ);
S = Si\I;

# a coefficients
a = zeros(q,1)*0im;
for i = 1:q
    a[i] = sum(S[i,:]);
end

# b coeeficitents
b = zeros(q,2*q)*0im;
SMinv = S*mi;
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
            c[j,m] = c[j,m] + li_hat[i+1](τ_t[m])*Si[i,j];
        end
    end
end
c0 = zeros(2*q,1);
for m = 1:2*q
   c0[m] = li_hat[1](τ_t[m]); 
end


# compute LHS and LU-decomposition
tid = time();
println("start")

M1 = φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ
save("LU_This.jld","M1",M1);

M2 = φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ
M3 = φᵢφⱼ + 1im*ϵ*k*γ[3]*∇φᵢ∇φⱼ + 1im*k*γ[3]*Vφᵢφⱼ 

save("LU_This.jld","M1",M1,"M2",M2,"M3");
println("DOOOOOOOONE")
LU_Matrix_1 = lu(φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ ) ;
LU_Matrix_2 = lu(φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ );
LU_Matrix_3 = lu(φᵢφⱼ + 1im*ϵ*k*γ[3]*∇φᵢ∇φⱼ + 1im*k*γ[3]*Vφᵢφⱼ );

# time integration
# set intial starting point for fixed point iteration
ξ_old_1 = zeros(SpaceDim_C)*0im;
ξ_old_2 = zeros(SpaceDim_C)*0im;
ξ_old_3 = zeros(SpaceDim_C)*0im;

U = zeros(SpaceDim_C)*0im;
copyto!(U,U0)


tid = time();




for n = 1:Nt
    
    # fixed-point itreration
    Δ = 1;
    it = 0;

    copyto!(ξ_old_1,U);
    copyto!(ξ_old_2,U);
    copyto!(ξ_old_3,U);
     tid = time();   
    while Δ > TOL && it<max_it
  
#  tid1 = time();      
        U_interim = c0[1]*U+c[1,1]*ξ_old_1+c[2,1]*ξ_old_2+c[3,1]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f1 = Assemble.u²uvᵢ(ω,U_interim,ρ)
 
        
        U_interim = c0[2]*U+c[1,2]*ξ_old_1+c[2,2]*ξ_old_2+c[3,2]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f2 = Assemble.u²uvᵢ(ω,U_interim,ρ)
        
        U_interim = c0[3]*U+c[1,3]*ξ_old_1+c[2,3]*ξ_old_2+c[3,3]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f3 = Assemble.u²uvᵢ(ω,U_interim,ρ)
        
        U_interim = c0[4]*U+c[1,4]*ξ_old_1+c[2,4]*ξ_old_2+c[3,4]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f4 = Assemble.u²uvᵢ(ω,U_interim,ρ)
  
      
        U_interim = c0[5]*U+c[1,5]*ξ_old_1+c[2,5]*ξ_old_2+c[3,5]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f5 = Assemble.u²uvᵢ(ω,U_interim,ρ)
  
        U_interim = c0[6]*U+c[1,6]*ξ_old_1+c[2,6]*ξ_old_2+c[3,6]*ξ_old_3;
        ρφᵢ = Assemble.ρvᵢ(U_interim,ω) 
        ρ = φᵢφⱼ_lu\ρφᵢ;
        f6 = Assemble.u²uvᵢ(ω,U_interim,ρ)
 # println("TIME ASSEMBLY ", time()-tid1)
  #tid2 = time();
        # compute right hand side
        F1	= a[1]*φᵢφⱼ*U - 1im*k*β*( b[1,1]*f1 + b[1,2]*f2 + b[1,3]*f3 +b[1,4]*f4+b[1,5]*f5+b[1,6]*f6 );                                       
        F2	= a[2]*φᵢφⱼ*U - 1im*k*β*( b[2,1]*f1 + b[2,2]*f2 + b[2,3]*f3 +b[2,4]*f4+b[2,5]*f5+b[2,6]*f6 );
        F3 	= a[3]*φᵢφⱼ*U - 1im*k*β*( b[3,1]*f1 + b[3,2]*f2 + b[3,3]*f3 +b[3,4]*f4+b[3,5]*f5+b[3,6]*f6 );
 	
 	ξ_new_1 = LU_Matrix_1\F1;
        ξ_new_2 = LU_Matrix_2\F2;
        ξ_new_3 = LU_Matrix_3\F3;
        # update iteration paramter
 #  println("TIME SOLVE ", time()-tid2);     
          
        ΔU_1 = ξ_old_1 - ξ_new_1; ΔU_2 = ξ_old_2 - ξ_new_2; ΔU_3 = ξ_old_3 - ξ_new_3
        Δ = sqrt(max( real(dot(ΔU_1,φᵢφⱼ*ΔU_1)),real(dot(ΔU_2,φᵢφⱼ*ΔU_2)),real(dot(ΔU_3,φᵢφⱼ*ΔU_3))) ) 
   
        copyto!(ξ_old_1,ξ_new_1);
        copyto!(ξ_old_2,ξ_new_2);
        copyto!(ξ_old_3,ξ_new_3);
   
        
        it+=1;

    end
       println(time()-tid, " it ", it)
 #  println(it)
  
    # set result and calculate U(n+1)
   # ξ_transform = change_basis(Si,ξ_new,Nx,q);
    ξ_transform_1 = Si[1,1]*ξ_old_1+Si[1,2]*ξ_old_2+Si[1,3]*ξ_old_3;
    ξ_transform_2 = Si[2,1]*ξ_old_1+Si[2,2]*ξ_old_2+Si[2,3]*ξ_old_3;
    ξ_transform_3 = Si[3,1]*ξ_old_1+Si[3,2]*ξ_old_2+Si[3,3]*ξ_old_3;
    
       
    

    U = li_hat[1](1)*U + li_hat[2](1)*ξ_transform_1+ li_hat[3](1)*ξ_transform_2+li_hat[4](1)*ξ_transform_3;
    
          if(mod(n,save_it)==0)
		
		ρφᵢ =Assemble.ρvᵢ(U,ω); 
 	#	println(sum(ρ_U))
	 	ρ = φᵢφⱼ_lu\ρφᵢ;
			 	
	 	E_nl = β/2*dot(ρ,φᵢφⱼ*ρ);
	 	#println(E_nl);
		Energy = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+E_nl;

		println("ENERGY"," ", (round((Energy)*10^7)/10^7).re,  " time ", n*k);
		
		save(save_here*"sol_SLOD_"*string(n)*".jld","U",U)
		 
	end

   
end
return U
end



function CG_q3(U0,T,Nt,ϕ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,φᵢφⱼ_lu,TOL,max_it,β,ϵ,save_it,save_here,Quad)
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


tid = time();



# compute LHS and LU-decomposition
LU_Matrix_1 = lu(φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ ) ;
LU_Matrix_2 = lu(φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ );
LU_Matrix_3 = lu(φᵢφⱼ + 1im*ϵ*k*γ[3]*∇φᵢ∇φⱼ + 1im*k*γ[3]*Vφᵢφⱼ );

# time integration
# set intial starting point for fixed point iteration

ξⁱ_1 = zeros(size(U0))*0im;
ξⁱ_2 = zeros(size(U0))*0im;
ξⁱ_3 = zeros(size(U0))*0im;

U = zeros(size(U0))*0im;
copyto!(U,U0)






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
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ;
        f₁ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
 
        
        Û = c0[2]*U+c[1,2]*ξⁱ_1+c[2,2]*ξⁱ_2+c[3,2]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₂ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
        
        Û = c0[3]*U+c[1,3]*ξⁱ_1+c[2,3]*ξⁱ_2+c[3,3]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₃ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
        
        Û = c0[4]*U+c[1,4]*ξⁱ_1+c[2,4]*ξⁱ_2+c[3,4]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
	ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₄ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
  
      
        Û = c0[5]*U+c[1,5]*ξⁱ_1+c[2,5]*ξⁱ_2+c[3,5]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₅ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
  
        Û = c0[6]*U+c[1,6]*ξⁱ_1+c[2,6]*ξⁱ_2+c[3,6]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ)) 
        ρₕ = ϕ'*ρ; 
        f₆ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
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
    
       
    

    U = li_hat[1](1)*U + li_hat[2](1)*ξ_transform_1+ li_hat[3](1)*ξ_transform_2+li_hat[4](1)*ξ_transform_3;
    
          if(mod(n,save_it)==0)
		
		Uₕ = ϕ'*U;
		ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Uₕ))
			 	
	 	E_nl = β/2*dot(ρ,φᵢφⱼ*ρ);
	 	#println(E_nl);
		Energy = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+E_nl;

		println("ENERGY"," ", (round((Energy)*10^7)/10^7).re,  " time ", n*k);
		
		save(save_here*"sol_SLOD_"*string(n)*".jld","U",U)
		 
	end

   
end
return U
end




function CG_q3_IterativeSolvers(U0,T,Nt,ϕ,∇φᵢ∇φⱼ,φᵢφⱼ,Vφᵢφⱼ,φᵢφ_lu,TOL,max_it,β,ϵ,save_it,save_here,Quad)
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


tid = time();

M1 = φᵢφⱼ + 1im*ϵ*k*γ[1]*∇φᵢ∇φⱼ + 1im*k*γ[1]*Vφᵢφⱼ
M2 = φᵢφⱼ + 1im*ϵ*k*γ[2]*∇φᵢ∇φⱼ + 1im*k*γ[2]*Vφᵢφⱼ
M3 = φᵢφⱼ + 1im*ϵ*k*γ[3]*∇φᵢ∇φⱼ + 1im*k*γ[3]*Vφᵢφⱼ 



# time integration
# set intial starting point for fixed point iteration

ξⁱ_1 = zeros(size(U0))*0im;
ξⁱ_2 = zeros(size(U0))*0im;
ξⁱ_3 = zeros(size(U0))*0im;

U = zeros(size(U0))*0im;
ρ = zeros(size(U0));
ξⁱ⁺¹_1 = zeros(size(U0))*0im;
ξⁱ⁺¹_2 = zeros(size(U0))*0im;
ξⁱ⁺¹_3 = zeros(size(U0))*0im;
copyto!(U,U0)

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
        #cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8 )
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ;
        f₁ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
 
        
        Û = c0[2]*U+c[1,2]*ξⁱ_1+c[2,2]*ξⁱ_2+c[3,2]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        #cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8)
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ; 
        f₂ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
        
        Û = c0[3]*U+c[1,3]*ξⁱ_1+c[2,3]*ξⁱ_2+c[3,3]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        #cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8)
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ; 
        f₃ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
        
        Û = c0[4]*U+c[1,4]*ξⁱ_1+c[2,4]*ξⁱ_2+c[3,4]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
	#cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8)
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ; 
        f₄ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
  
      
        Û = c0[5]*U+c[1,5]*ξⁱ_1+c[2,5]*ξⁱ_2+c[3,5]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        #cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8)
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ; 
        f₅ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
  
        Û = c0[6]*U+c[1,6]*ξⁱ_1+c[2,6]*ξⁱ_2+c[3,6]*ξⁱ_3;
        Ûₕ = ϕ'*Û;
        #cg!(ρ,φᵢφⱼ,ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ),reltol =1e-8)
        ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Ûₕ))
        ρₕ = ϕ'*ρ; 
        f₆ = ϕ*Assemble.ρuvᵢ(mesh,Quad,Ûₕ,ρₕ);
 # println("TIME ASSEMBLY ", time()-tid1)
  #tid2 = time();
        # compute right hand side
        F1	= a[1]*φᵢφⱼ*U - 1im*k*β*( b[1,1]*f₁ + b[1,2]*f₂ + b[1,3]*f₃ +b[1,4]*f₄+b[1,5]*f₅+b[1,6]*f₆ );                                       
        F2	= a[2]*φᵢφⱼ*U - 1im*k*β*( b[2,1]*f₁ + b[2,2]*f₂ + b[2,3]*f₃ +b[2,4]*f₄+b[2,5]*f₅+b[2,6]*f₆ );
        F3 	= a[3]*φᵢφⱼ*U - 1im*k*β*( b[3,1]*f₁ + b[3,2]*f₂ + b[3,3]*f₃ +b[3,4]*f₄+b[3,5]*f₅+b[3,6]*f₆ );
 	
#@time 	gmres!(ξⁱ⁺¹_1,M1,F1);#gauss_seidel!(ξⁱ⁺¹_1,M1,F1,maxiter=10);#gmres!(ξⁱ⁺¹_1,M1,F1,reltol =1e-8); #LU_Matrix_1\F1;
println("------ gmres ------")
@time	ξⁱ⁺¹_1=gmres(M1,F1);
 	#gmres!(ξⁱ⁺¹_2,M2,F2);#gauss_seidel!(ξⁱ⁺¹_2,M2,F2,maxiter=10);#gmres!(ξⁱ⁺¹_2,M2,F2,reltol =1e-8);
@time	ξⁱ⁺¹_2=gmres(M2,F2);
 	#gmres!(ξⁱ⁺¹_3,M3,F3);#gauss_seidel!(ξⁱ⁺¹_3,M3,F3,maxiter=10);#gmres!(ξⁱ⁺¹_3,M3,F3,reltol =1e-8);
@time 	ξⁱ⁺¹_3=gmres(M3,F3);
 	
 	#ξⁱ⁺¹_1 = LU_Matrix_1\F1;
        #ξⁱ⁺¹_2 = LU_Matrix_2\F2;
        #ξⁱ⁺¹_3 = LU_Matrix_3\F3;
        
        # update iteration paramter
 #  println("TIME SOLVE ", time()-tid2);     
          
        ΔU_1 = ξⁱ_1 - ξⁱ⁺¹_1; ΔU_2 = ξⁱ_2 - ξⁱ⁺¹_2; ΔU_3 = ξⁱ_3 - ξⁱ⁺¹_3
        Δ = sqrt(max( real(dot(ΔU_1,φᵢφⱼ*ΔU_1)),real(dot(ΔU_2,φᵢφⱼ*ΔU_2)),real(dot(ΔU_3,φᵢφⱼ*ΔU_3))) ) 
   	println(Δ);
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
    
       
    

    U = li_hat[1](1)*U + li_hat[2](1)*ξ_transform_1+ li_hat[3](1)*ξ_transform_2+li_hat[4](1)*ξ_transform_3;
    
          if(mod(n,save_it)==0)
		
		Uₕ = ϕ'*U;
		ρ = φᵢφⱼ_lu\(ϕ*Assemble.u²vᵢ(mesh,Quad,Uₕ))
			 	
	 	E_nl = β/2*dot(ρ,φᵢφⱼ*ρ);
	 	#println(E_nl);
		Energy = U'*(ϵ*(∇φᵢ∇φⱼ*U)+(Vφᵢφⱼ*U))+E_nl;

		println("ENERGY"," ", (round((Energy)*10^7)/10^7).re,  " time ", n*k);
		
		save(save_here*"sol_SLOD_"*string(n)*".jld","U",U)
		 
	end

   
end
return U
end

