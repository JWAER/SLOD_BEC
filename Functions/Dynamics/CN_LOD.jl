function  CN_LOD(U0,k,Nt,A_LOD,M_LOD,MV_LOD,W,W_prim,TOL,max_it,β,ϵ,save_it,Vt,dofs_f,VLOD);



U = zeros(SpaceDim_C)*0im;
Unext = zeros(SpaceDim_C)*0im;
copyto!(U,U0);

		ρ = Assemble_Rho(U,W) 
      		ρ_LOD=Mlu\ρ;
      
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+β/2*ρ_LOD'*(M_LOD*ρ_LOD) 
		println("ENERGY"," ", round(real(Energy)*10^7)/10^7);

TRI = sparse(LowerTriangular(M_LOD));

PDE_Matrix_LU = lu(M_LOD+1im*k/2*(ϵ*A_LOD+MV_LOD));

for idx_t = 1:Nt

		Δ = 1;
		it = 0;


		tm = (idx_t-1/2)*k;
		Vt_rhs = Compute_RHS(Vt,tm,Quad_6,Mesh,dofs_f) 
		Vt_LOD = Mlu\(VLOD*Vt_rhs);
	
		MVt_LOD = Assemble_M_NL(Vt_LOD,W,W_prim,TRI)
		MVt_LOD += MVt_LOD'; 
		for d = 1:SpaceDim_C; MVt_LOD[d,d]*=1/2;end
	
		PDE_Matrix_LOD = M_LOD+1im*k/2*(ϵ*A_LOD+MVt_LOD+MV_LOD);

		
		copyto!(Unext,U);
		CONST = (M_LOD-1im*k/2*(ϵ*A_LOD+MV_LOD+MVt_LOD))*U;
	
		ρ0 = Assemble_Rho(U,W);

#save("hm.jld", "Pde",PDE_Matrix_LOD,"MVt_LOD",MVt_LOD); println("Done")
#pdeLU = lu(PDE_Matrix_LOD);
	 while( Δ > 10.0^-9 && it < 15)
	 	
	 	Um = (U+Unext)/2;
	
		ρ1 = Assemble_Rho(Unext,W) 		 	
		ρm = (ρ0+ρ1)/2
	 	
	 	ρm_LOD = Mlu\ρm	

	 	NL = Assemble_NL(W,Um,ρm_LOD)
        

		rhs_LOD = CONST-β*1im*k*NL
	
		#Unext = PDE_Matrix_LU\f_VLOD;#
	
		#Unext = cg(PDE_Matrix_LOD,rhs_LOD,reltol=10^-14);
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD)
		#Unext=PDE_Matrix_LOD\rhs_LOD
		
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU);
		#println(norm(PDE_Matrix_LOD*Unext-rhs_LOD))
		
#		Unext = PDE_Matrix_LOD\rhs_LOD
#@time		hm = PDE_Matrix_LOD\rhs_LOD
@time		gmres!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#		cg!(Unext,	PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#@time		Unext = pdeLU\rhs_LOD;

#		copyto!(Unext,hm);
#		println("2 ", norm(PDE_Matrix_LOD*Unext1-rhs_LOD)," ",norm(PDE_Matrix_LOD*Unext-rhs_LOD), " ",norm(Unext-Unext1));
		#copyto!(Unext,Unext1);
		
		
		ΔU = (2*Um-U-Unext);
		Δ = sqrt(real(dot(ΔU,M_LOD*ΔU)))   
		       it+=1
		       println(Δ)



	end
	println(it)
	 copyto!(U,Unext);

		if(mod(idx_t,save_it)==0)
		
		ρ_interim =Assemble_Rho(U,W); 
	 	ρ_LOD_interim=Mlu\ρ_interim;
			 	
	 	E_nl = β/2*dot(ρ_LOD_interim,M_LOD*ρ_LOD_interim);
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+E_nl;

		println("ENERGY"," ", round((Energy)*10^7)/10^7,  " time ", idx_t*k, " mass ", dot(U,M_LOD,U));
		
		save("./Results/Dynamics/sol_LOD_"*string(idx_t)*".jld","U",U)
		 
		end


	end

return U
end





















function  CN_LOD_Imp(U0,k,Nt,A_LOD,M_LOD,MV_LOD,W,W_prim,TOL,max_it,β,ϵ,save_it,Vt,dofs_f,VLOD);



U = zeros(SpaceDim_C)*0im;
Unext = zeros(SpaceDim_C)*0im;
copyto!(U,U0);

		ρ = Assemble_Rho(U,W) 
      		ρ_LOD=Mlu\ρ;
      
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+β/2*ρ_LOD'*(M_LOD*ρ_LOD) 
		println("ENERGY"," ", round(real(Energy)*10^7)/10^7);

TRI = sparse(LowerTriangular(M_LOD));

PDE_Matrix_LU = lu(M_LOD+1im*k/2*(ϵ*A_LOD+MV_LOD));

for idx_t = 1:Nt

		Δ = 1;
		it = 0;


		tm = (idx_t-1/2)*k;
		Vt_rhs = Compute_RHS(Vt,tm,Quad_6,Mesh,dofs_f) 
		Vt_LOD = Mlu\(VLOD*Vt_rhs);
	
		MVt_LOD = Assemble_M_NL(Vt_LOD,W,W_prim,TRI)
		MVt_LOD += MVt_LOD'; 
		for d = 1:SpaceDim_C; MVt_LOD[d,d]*=1/2;end

		ρ0 = Assemble_Rho(U,W);
		ρ0_LOD = Mlu\ρ0

		Mρ_LOD = Assemble_M_NL(ρ0_LOD,W,W_prim,TRI)
		Mρ_LOD += Mρ_LOD'; 
		for d = 1:SpaceDim_C; Mρ_LOD[d,d]*=1/2;end
		PDE_Matrix_LOD = M_LOD+1im*k/2*(ϵ*A_LOD+MVt_LOD+MV_LOD+Mρ_LOD/2);

		
		copyto!(Unext,U);
		CONST = (M_LOD-1im*k/2*(ϵ*A_LOD+MV_LOD+MVt_LOD+Mρ_LOD/2))*U;
	


#save("hm.jld", "Pde",PDE_Matrix_LOD,"MVt_LOD",MVt_LOD); println("Done")
tid = time();

	 while( Δ > 10.0^-9 && it < 15)
	 	
	 	Um = (U+Unext)/2;
	
		ρ1 = Assemble_Rho(Unext,W) 		 	
		
	 	
	 	ρ1_LOD = Mlu\ρ1	

	 	NL = Assemble_NL(W,Um,ρ1_LOD)
        

		rhs_LOD = CONST-β*1im*k*NL/2
	
		#Unext = PDE_Matrix_LU\f_VLOD;#
	
		#Unext = cg(PDE_Matrix_LOD,rhs_LOD,reltol=10^-14);
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD)
		#Unext=PDE_Matrix_LOD\rhs_LOD
		
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU);
		#println(norm(PDE_Matrix_LOD*Unext-rhs_LOD))
		
#		Unext = PDE_Matrix_LOD\rhs_LOD
#@time		hm = PDE_Matrix_LOD\rhs_LOD
@time		gmres!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#		cg!(Unext,	PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#@time		Unext = pdeLU\rhs_LOD;

#		copyto!(Unext,hm);
#		println("2 ", norm(PDE_Matrix_LOD*Unext1-rhs_LOD)," ",norm(PDE_Matrix_LOD*Unext-rhs_LOD), " ",norm(Unext-Unext1));
		#copyto!(Unext,Unext1);
		
		
		ΔU = (2*Um-U-Unext);
		Δ = sqrt(real(dot(ΔU,M_LOD*ΔU)))   
		println(Δ)
		       it+=1


	end
	println(time() -tid, " ",it)
	 copyto!(U,Unext);

		if(mod(idx_t,save_it)==0)
		
		ρ_interim =Assemble_Rho(U,W); 
	 	ρ_LOD_interim=Mlu\ρ_interim;
			 	
	 	E_nl = β/2*dot(ρ_LOD_interim,M_LOD*ρ_LOD_interim);
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+E_nl;

		println("ENERGY"," ", round((Energy)*10^7)/10^7,  " time ", idx_t*k, " mass ", dot(U,M_LOD,U));
		
		save("./Results/Dynamics/sol_LOD_"*string(idx_t)*".jld","U",U)
		 
		end


	end

return U
end




function  CN_LOD_exp(U0,k,Nt,A_LOD,M_LOD,MV_LOD,W,W_prim,TOL,max_it,β,ϵ,save_it,Vt,dofs_f,VLOD);



U = zeros(SpaceDim_C)*0im;
Unext = zeros(SpaceDim_C)*0im;
copyto!(U,U0);

		ρ = Assemble_Rho(U,W) 
      		ρ_LOD=Mlu\ρ;
      
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+β/2*ρ_LOD'*(M_LOD*ρ_LOD) 
		println("ENERGY"," ", round(real(Energy)*10^7)/10^7);

TRI = sparse(LowerTriangular(M_LOD));

PDE_Matrix_LU = lu(M_LOD+1im*k/2*(ϵ*A_LOD+MV_LOD));

for idx_t = 1:Nt
tid = time()
		Δ = 1;
		it = 0;
TID = time()

		tm = (idx_t-1/2)*k;
		Vt_rhs = Compute_RHS(Vt,tm,Quad_6,Mesh,dofs_f) 
		Vt_LOD = Mlu\(VLOD*Vt_rhs);
	
		MVt_LOD = Assemble_M_NL(Vt_LOD,W,W_prim,TRI)
		MVt_LOD += MVt_LOD'; 
		for d = 1:SpaceDim_C; MVt_LOD[d,d]*=1/2;end
	
		PDE_Matrix_LOD = M_LOD+1im*k/2*(ϵ*A_LOD+MVt_LOD+MV_LOD);

		
		copyto!(Unext,U);
		CONST = M_LOD*U-1im*k/2*(ϵ*A_LOD*U+MV_LOD*U);
	
		ρ0 = Assemble_Rho(U,W);
println("outside loop ", time()-TID);
#save("hm.jld", "Pde",PDE_Matrix_LOD,"MVt_LOD",MVt_LOD); println("Done")
#pdeLU = lu(PDE_Matrix_LOD);
TID = time();
	 while( Δ > 10.0^-9 && it < 15)
	 	
	 	Um = (U+Unext)/2;
	
		ρ1 = Assemble_Rho(Unext,W) 		 	
		ρm = (ρ0+ρ1)/2
	 	
	 	ρm_LOD = Mlu\ρm	
println("nl")
	 @time	NL = Assemble_NL(W,Um,ρm_LOD)
        

		rhs_LOD = CONST-β*1im*k*NL-1im*k*MVt_LOD*Um
	
		#Unext = PDE_Matrix_LU\f_VLOD;#
	
		#Unext = cg(PDE_Matrix_LOD,rhs_LOD,reltol=10^-14);
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD)
		#Unext=PDE_Matrix_LOD\rhs_LOD
		
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU);
		#println(norm(PDE_Matrix_LOD*Unext-rhs_LOD))
		
#		Unext = PDE_Matrix_LOD\rhs_LOD
#@time		hm = PDE_Matrix_LOD\rhs_LOD
#@time		gmres!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#		cg!(Unext,	PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
@time		Unext = PDE_Matrix_LU\rhs_LOD;

#		copyto!(Unext,hm);
#		println("2 ", norm(PDE_Matrix_LOD*Unext1-rhs_LOD)," ",norm(PDE_Matrix_LOD*Unext-rhs_LOD), " ",norm(Unext-Unext1));
		#copyto!(Unext,Unext1);
		
		
		ΔU = (2*Um-U-Unext);
		Δ = sqrt(real(dot(ΔU,M_LOD*ΔU)))   
#		println(Δ)
		       it+=1


	end
#	println(it," ", Δ, " " , time()-TID)
	 copyto!(U,Unext);

		if(mod(idx_t,save_it)==0)
		
		ρ_interim =Assemble_Rho(U,W); 
	 	ρ_LOD_interim=Mlu\ρ_interim;
			 	
	 	E_nl = β/2*dot(ρ_LOD_interim,M_LOD*ρ_LOD_interim);
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+E_nl;

		println("ENERGY"," ", round((Energy)*10^7)/10^7,  " time ", idx_t*k, " mass ", dot(U,M_LOD,U));
		
		save("./Results/Dynamics/sol_LOD_"*string(idx_t)*".jld","U",U)
		 
		end


	end

return U
end





function  CN_LOD_exp_NO_DYN(Mesh,U0,k,Nt,A_LOD,M_LOD,MV_LOD,W_h,W_h_prim,TOL,max_it,β,ϵ,save_it,dofs_f,VLOD,Mlu);

SpaceDim_C = size(A_LOD,2);

U = zeros(SpaceDim_C)*0im;
Unext = zeros(SpaceDim_C)*0im;
copyto!(U,U0);

		U_h = VLOD'*U;
		ρ_LOD = Mlu\(VLOD*Assemble_Rho(U_h,W_h))
      
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+β/2*ρ_LOD'*(M_LOD*ρ_LOD) 
		println("ENERGY"," ", round(real(Energy)*10^7)/10^7);

#TRI = sparse(LowerTriangular(M_LOD));

PDE_Matrix_LU = lu(M_LOD+1im*k/2*(ϵ*A_LOD+MV_LOD));
PDE_Matrix_LOD = M_LOD+1im*k/2*(ϵ*A_LOD+MV_LOD);

tiid = time();

for idx_t = 1:Nt
tt = time();

#println("NEW STEP-----------------------------------")
		Δ = 1;
		it = 0;
#=
		tm = (idx_t-1/2)*k;
		Vt_rhs = Compute_RHS(Vt,tm,Quad_6,Mesh,dofs_f) 
		Vt_LOD = Mlu\(VLOD*Vt_rhs);
	
		MVt_LOD = Assemble_M_NL(Vt_LOD,W,W_prim,TRI)
		MVt_LOD += MVt_LOD'; 
		for d = 1:SpaceDim_C; MVt_LOD[d,d]*=1/2;end
=#	
	
		
		copyto!(Unext,U);
		CONST = M_LOD*U-1im*k/2*(ϵ*A_LOD*U+MV_LOD*U);
		
		U_h = VLOD'*U;
		ρ0 = VLOD*Assemble_Rho(U_h,W_h);
		
#println("outside loop ", time()-TID);
#save("hm.jld", "Pde",PDE_Matrix_LOD,"MVt_LOD",MVt_LOD); println("Done")
#pdeLU = lu(PDE_Matrix_LOD);
	 while( Δ > TOL && it < 15)
	 	Um = (U+Unext)/2;
	 	Um_h = VLOD'*Um
	 	
	 	Unext_h = VLOD'*Unext;
		ρ1 = VLOD*Assemble_Rho(Unext_h,W_h) 		 	
		ρm = (ρ0+ρ1)/2
	 	
	 	ρm_LOD = Mlu\ρm
	 	ρm_h = VLOD'*ρm_LOD	
	 	NL = VLOD*Assemble_NL(W_h,Um_h,ρm_h)
        

		rhs_LOD = CONST-β*1im*k*NL
	
		#Unext = PDE_Matrix_LU\f_VLOD;#
	
		#Unext = cg(PDE_Matrix_LOD,rhs_LOD,reltol=10^-14);
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD)
		#Unext=PDE_Matrix_LOD\rhs_LOD
		
		#cg!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU);
		#println(norm(PDE_Matrix_LOD*Unext-rhs_LOD))
		
#		Unext = PDE_Matrix_LOD\rhs_LOD
#@time		hm = PDE_Matrix_LOD\rhs_LOD
#@time		gmres!(Unext,PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
#		cg!(Unext,	PDE_Matrix_LOD,rhs_LOD,Pl=PDE_Matrix_LU)
		Unext = PDE_Matrix_LU\rhs_LOD;
		
		ΔU = (2*Um-U-Unext);
		Δ = sqrt(real(dot(ΔU,M_LOD*ΔU)))   
#		println(Δ)
		       it+=1


	end
	println(time()-tt, " ", it)
#	println(it," ", Δ, " " , time()-TID)
	 copyto!(U,Unext);

		if(idx_t==Nt)
		 
		U_h = VLOD'*U;
	 	ρ_LOD_interim=Mlu\(VLOD*Assemble_Rho(U_h,W_h););
			 	
	 	E_nl = β/2*dot(ρ_LOD_interim,M_LOD*ρ_LOD_interim);
		Energy = U'*(ϵ*(A_LOD*U)+(MV_LOD*U))+E_nl;

		println("ENERGY"," ", round((Energy)*10^7)/10^7,  " time ", idx_t*k, " mass ", dot(U,M_LOD,U));
		
#		save("./Results/Dynamics/sol_LOD_"*string(idx_t)*".jld","U",U)
		 
		end


	end
println("online time ", time()-tiid)
return U
end



