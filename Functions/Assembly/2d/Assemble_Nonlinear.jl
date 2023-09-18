function fₕvᵢvⱼ(mesh,Quad,f,sparsity) #f must be real, returns real valued matrix!

	M = similar(sparsity); fill!(M.nzval,0.);
	
	Dim = size(mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);
	M_loc = zeros(10,10);
	
        X_gp = Quad["X"];
        Y_gp = Quad["Y"];
        ω_gp = Quad["ω"];

        N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
      ω_gp = reshape(ω_gp,1,N_gp);
	#-----------------------------------------------------------------------------------------   

	quad_vals = zeros(10);

	W_loc = zeros(10,10,10);

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	map_f = zeros(Int,size(mesh.p,2));
	for i = 1:length(map_f); if(!isempty(searchsorted(dofs_f,i)));map_f[i] = searchsorted(dofs_f,i)[1]; end; end

	loc2glob_dirichlet = map_f[mesh.loc2glob];

	f_loc = zeros(10);
	
	loc2glob = zeros(Int,10);
	simplex = zeros(Int,3);
	Vertices = zeros(2,3);

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)
	fill!(f_loc,0);
	fill!(M_loc,0.);
	
	simplex .= mesh.t_h[:,n_t];
	Vertices .=  mesh.p[:,simplex];
	loc2glob .=loc2glob_dirichlet[:,n_t]
	
	
      	detJ = det(Vertices[:,2:3].-Vertices[:,1]);
      	
	for i = 1:10; if(loc2glob[i]!=0); f_loc[i] = f[loc2glob[i]]; end end ;

	
	
		if(n_t==1)	#computes reference tensor
			for gp in 1:N_gp
				fill!(quad_vals,0);
	
				idx = 1;
			
			        #println(len)

			        for i = 1:10; quad_vals[i] = vₕ[i](X_gp[gp],Y_gp[gp])	end	
				for i = 1:10
					for j = 1:10
						for k = 1:10
							W_loc[i,j,k] += quad_vals[i]*quad_vals[j]*quad_vals[k]*ω_gp[gp];
						end
					end
				end
			end
		end
	#tid_loc = time();
		for i = 1:10
			for j = i:10
				c = 0;
				for k = 1:10
					c+= W_loc[k,j,i]*f_loc[k] # result is symmetric in i,j
				end
				M_loc[j,i] = c*detJ; c*=0; 
				
			end
		end
		M_loc+=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
	#	println(M_loc)
	#println("loc ", time()-tid_loc)
	#tid_add = time();
	#add to global tensor using local tensor which is constant!
		for j = 1:10; j_glob = loc2glob[j]; if(j_glob!=0);
			idx = M.colptr[j_glob]:M.colptr[j_glob+1]-1;
			row = M.rowval[idx];
			for i=1:10; i_glob = loc2glob[i]; if(i_glob!=0);
		#		println(i_glob, " ", j_glob, " ", M_loc[i,j]);
				m = searchsorted(row,i_glob)[1];
				M.nzval[idx[m]]+=M_loc[j,i];
		#		M[i_glob,j_glob]+=M_loc[i,j];
		
			end end
		end end
	#println("add ", time()-tid_add);
	end
		return M 
			


end




function u²vᵢ(mesh,Quad,U)

	Rho = zeros(size(U));
	Rho_loc = zeros(10);
	
	N_gp = Int(Quad["N_gp"][1]);
	
        X_gp = Quad["X"];
        Y_gp = Quad["Y"];
        ω_gp = Quad["ω"];

        N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	ω_gp = reshape(ω_gp,1,N_gp);
	#-----------------------------------------------------------------------------------------   

	quad_vals = zeros(10);

	W_loc = zeros(10,10,10);

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	map_f = zeros(Int,size(mesh.p,2));
	for i = 1:length(map_f); if(!isempty(searchsorted(dofs_f,i)));map_f[i] = searchsorted(dofs_f,i)[1]; end; end

	loc2glob_dirichlet = map_f[mesh.loc2glob];

	ta = time();
	u_loc = zeros(typeof(U[1]),10);
	


	loc2glob = zeros(Int,10);
	simplex = zeros(Int,3);
	Vertices = zeros(2,3);

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)
	
	simplex .=  mesh.t_h[:,n_t]#mesh.t_h[[1,4,10],n_t];
	Vertices .=  mesh.p[:,simplex];
	loc2glob .=loc2glob_dirichlet[:,n_t]
	
	
      	detJ = det(Vertices[:,2:3].-Vertices[:,1]);#
    
	fill!(u_loc,0.)
	for i = 1:10; if(loc2glob[i]!=0); u_loc[i] = U[loc2glob[i]]; end end ;

	
	
		if(n_t==1)	#computes reference tensor
			for gp in 1:N_gp
				fill!(quad_vals,0);
	
				idx = 1;
			
			 
			        for i = 1:10; quad_vals[i] = vₕ[i](X_gp[gp],Y_gp[gp]) end	
				for i = 1:10
					for j = 1:10
						for k = 1:10
							W_loc[i,j,k] += quad_vals[i]*quad_vals[j]*quad_vals[k]*ω_gp[gp];
						end
					end
				end
			end
			
			for i = 1:10
				for j = 1:10
					W_loc[j,j,i]/=2;
				end
			end
		end
	
	#tid_loc = time();
		for i = 1:10
			for j = 1:10;
				ū = conj(u_loc[j]);
				for k = j:10
					Rho_loc[i]+= W_loc[k,j,i]*2*real(u_loc[k]*ū);
				end
			end
		end

	for i = 1:10; i_glob = loc2glob[i]; if(i_glob!=0);
		 Rho[i_glob] += Rho_loc[i]*detJ; Rho_loc[i] = 0;
		end end
	end
		return Rho 
			


end




function ρuvᵢ(mesh,Quad,U, ρ)

	NL = similar(U); fill!(NL,0);

	NL_loc = zeros(typeof(U[1]),10);
	
	u_loc = zeros(typeof(U[1]),10);
	ρ_loc = zeros(10);
	
	N_gp = Int(Quad["N_gp"][1]);
	
        X_gp = Quad["X"];
        Y_gp = Quad["Y"];
        ω_gp = Quad["ω"];

        N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	ω_gp = reshape(ω_gp,1,N_gp);
	#-----------------------------------------------------------------------------------------   

	quad_vals = zeros(10);

	W_loc = zeros(10,10,10);

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	map_f = zeros(Int,size(mesh.p,2));
	for i = 1:length(map_f); if(!isempty(searchsorted(dofs_f,i)));map_f[i] = searchsorted(dofs_f,i)[1]; end; end

	loc2glob_dirichlet = map_f[mesh.loc2glob];

	


	loc2glob = zeros(Int,10);
	simplex = zeros(Int,3);
	Vertices = zeros(2,3);

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)
	fill!(u_loc,0.)
	
	simplex .= mesh.t_h[:,n_t];
	Vertices .=  mesh.p[:,simplex];
	loc2glob .=loc2glob_dirichlet[:,n_t]
	
	
      	detJ = det(Vertices[:,2:3].-Vertices[:,1]);#
    

	for i = 1:10; if(loc2glob[i]!=0); u_loc[i] = U[loc2glob[i]];
	 ρ_loc[i] = ρ[loc2glob[i]] end end ;

	
	
		if(n_t==1)	#computes reference tensor
			for gp in 1:N_gp
				fill!(quad_vals,0);
	
				idx = 1;
			
			 
			        for i = 1:10; quad_vals[i] = vₕ[i](X_gp[gp],Y_gp[gp]) end	
				for i = 1:10
					for j = 1:10
						for k = 1:10
							W_loc[i,j,k] += quad_vals[i]*quad_vals[j]*quad_vals[k]*ω_gp[gp];
						end
					end
				end
			end
		end
	
	#tid_loc = time();
		for i = 1:10
			for j = 1:10;
				for k = 1:10
					NL_loc[i]+= W_loc[k,j,i]*u_loc[k]*ρ_loc[j];
				end
			end
		end

	for i = 1:10; i_glob = loc2glob[i]; if(i_glob!=0);
		 NL[i_glob] += NL_loc[i]*detJ; NL_loc[i] = 0;
		end end
	end
		return NL 
			


end


















