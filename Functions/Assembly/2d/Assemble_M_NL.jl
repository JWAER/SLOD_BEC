function fₕvᵢvⱼ(Mesh,Quad,f,M)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);
	M_loc = zeros(10,10);

	X_gp = [ 0.045189009784400
	   0.045189009784400
	   0.909621980431200
	   0.747512472733900
	   0.222063165537300
	   0.747512472733900
	   0.222063165537300
	   0.030424361728800
	   0.030424361728800
	   0.136991201264900
	   0.644718727763700
	   0.136991201264900
	   0.218290070971400
	   0.218290070971400
	   0.644718727763700
	   0.036960330433400
	   0.481519834783300
	   0.481519834783300
	   0.403603979817900
	   0.403603979817900
	   0.192792040364100];
	   
	Y_gp = [0.045189009784400
	   0.909621980431200
	   0.045189009784400
	   0.030424361728800
	   0.030424361728800
	   0.222063165537300
	   0.747512472733900
	   0.747512472733900
	   0.222063165537300
	   0.218290070971400
	   0.218290070971400
	   0.644718727763700
	   0.644718727763700
	   0.136991201264900
	   0.136991201264900
	   0.481519834783300
	   0.036960330433400
	   0.481519834783300
	   0.192792040364100
	   0.403603979817900
	   0.403603979817900];


	w_gp =  [0.051987142064600
	   0.051987142064600
	   0.051987142064600
	   0.070703410178400
	   0.070703410178400
	   0.070703410178400
	   0.070703410178400
	   0.070703410178400
	   0.070703410178400
	   0.090939076095200
	   0.090939076095200
	   0.090939076095200
	   0.090939076095200
	   0.090939076095200
	   0.090939076095200
	   0.103234405138000
	   0.103234405138000
	   0.103234405138000
	   0.188160146916700
	   0.188160146916700
	   0.188160146916700]/4;

	N_gp = length(w_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	w_gp = reshape(w_gp,1,N_gp);
	#-----------------------------------------------------------------------------------------   

	quad_vals = zeros(10);

	W_loc = zeros(10,10,10);

	dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
	map_f = zeros(Int,size(mesh.p,2));
	for i = 1:length(map_f); if(!isempty(searchsorted(dofs_f,i)));map_f[i] = searchsorted(dofs_f,i)[1]; end; end

	loc2glob_dirichlet = map_f[mesh.loc2glob];

	ta = time();
	
	f_loc = zeros(10);
	
#tid = time();

	loc2glob = zeros(Int,10);
	simplex = zeros(Int,3);
	Vertices = zeros(2,3);

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)
	fill!(f_loc,0);
	fill!(M_loc,0.);
	
		simplex .= mesh.t_h[:,n_t];
		Vertices .=  mesh.p[:,simplex];
		#J = hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1]);
		#Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp]; 
		loc2glob .=loc2glob_dirichlet[:,n_t]
		
		#detJ = det(J)
		
		detJ = (Vertices[1,2]-Vertices[1,1])*(Vertices[2,3]-Vertices[2,1])-(Vertices[1,3]-Vertices[1,1])*(Vertices[2,2]-Vertices[2,1]);
		

	for i = 1:10; if(loc2glob[i]!=0); f_loc[i] = f[loc2glob[i]]; end end ;

	
	
		if(n_t==1)	#computes reference tensor
			for gp in 1:N_gp
				fill!(quad_vals,0);
	
				idx = 1;
			
			        #println(len)

			        for i = 1:10; quad_vals[i] = vₕ[i](X_gp[gp],Y_gp[gp])	end	#calculate value in GP ! MESHIDX mustnt be zero!
				for i = 1:10
					for j = 1:10
						for k = 1:10
							W_loc[i,j,k] += quad_vals[i]*quad_vals[j]*quad_vals[k]*w_gp[gp];
						end
					end
				end
			end
		end
	tid_loc = time();
	
		for i = 1:10
			for j = 1:10
				c = 0;
				for k = 1:10
					c+= W_loc[i,j,k]*f_loc[k]
				end
				M_loc[i,j] = c*detJ; c*=0; 
				
			end
		end
	#	println(M_loc)
	#println("loc ", time()-tid_loc)
	#add to global tensor using local tensor which is constant!

	tid_add = time();
		for j = 1:10; j_glob = loc2glob[j]; if(j_glob!=0);
#			idx = M.colptr[j_glob]:M.colptr[j_glob+1]-1;
#			row = M.rowval[idx];
			for i=1:10; i_glob = loc2glob[i]; if(i_glob!=0);
				#m = searchsorted(row,i_glob)[1];
				#M.nzval[idx[m]]+=M_loc[i,j];
				M[i_glob,j_glob]+=M_loc[i,j];
		
			end end
		end end
	#println("add ", time()-tid_add);
	end
		return M 
			


end


