# W[i][j,k] beräknas för j,k>= i vilket betyder att W[i][j,k] för j,k<i, antag k<=j<i, sätt w[i][j,k] = w[k][j,i]

#For each coarse simplex find quadrature points, then find corresponding fine simplex and do as usual!
function  Compute_Wh(Mesh,ω,Quad)


	#----------------Quadrature rule to EXACTLY integrate a 9th degree polynomial-----------------
        X_gp = Quad["X"];
        Y_gp = Quad["Y"];
        Z_gp = Quad["Z"];
        ω_gp = Quad["ω"];

        N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
        Z_gp = reshape(Z_gp,1,N_gp);
	ω_gp = reshape(ω_gp,1,N_gp);
	#-----------------------------------------------------------------------------------------   

	quad_vals = zeros(20);

	W_loc = zeros(20,20,20);

map_f = zeros(Int,size(Mesh.p,2));
for i = 1:length(map_f); if(!isempty(searchsorted(Mesh.dofs,i)));map_f[i] = searchsorted(Mesh.dofs,i)[1]; end; end

loc2glob_dirichlet = map_f[Mesh.t];

	ta = time();
	
#tid = time();

#W_loc = load("W_loc_P3.jld")["W_loc"];

	for n_t = 1:size(Mesh.t,2)#size(t_C,2)
	        Simplex = Mesh.t[[1,4,10,20],n_t];
		Vertices =  Mesh.p[:,Simplex];
                J = Vertices[:,2:end].-Vertices[:,1];
		Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp;Z_gp]; 
		
		loc2glob =loc2glob_dirichlet[:,n_t] #dont do setdiff(...,0), will destroy correspondance with W_loc
	
	if(n_t==1)	
		for gp in 1:N_gp
			fill!(quad_vals,0);

			idx = 1;
			
		        #println(len)

                        for i = 1:20; quad_vals[i] = vₕ[i]([X_gp[gp],Y_gp[gp],Z_gp[gp]])	end	#calculate value in GP ! MESHIDX mustnt be zero!
			for i = 1:20
				for j = 1:20
					for k = 1:20
						W_loc[i,j,k] += quad_vals[i]*quad_vals[j]*quad_vals[k]*ω_gp[gp];
					end
				end
			end
		end
	end
	#add to global tensor using local tensor which is constant!

	detJ = det(J)

	Quad_Rule(W_loc,detJ,loc2glob,ω) 	
	
		
	end 


end


function Quad_Rule(W_loc,detJ,loc2glob,ω)
#sort loc2glob? or use if-statements

	for i = 1:20
		i_glob = loc2glob[i];
	if(i_glob!=0)
		sub_i = ω.Iptr[i_glob]: ω.Iptr[i_glob+1]-1;

		PreAlloc_IdxJ = ω.J[sub_i];

		for j = 1:20
			j_glob = loc2glob[j]
		if(i_glob<=j_glob);	
			sub_ij= (sub_i[1]-1).+searchsorted(PreAlloc_IdxJ, j_glob)

		if(!isempty(sub_ij))	
	
			PreAlloc_IdxK = ω.K[sub_ij]
			for k = 1:20
			k_glob = loc2glob[k];
			if(j_glob<=k_glob)
				sub_ijk=(sub_ij[1]-1).+searchsorted(PreAlloc_IdxK, k_glob)
				if(!isempty(sub_ijk)) ω.Val[sub_ijk[1]]+=W_loc[i,j,k]*detJ  end
			end
			end
		end
		end
		end

	end
	end
end











#---------------------------------------------#
#=
No 21 is in correct (2,3,145) = 0.00019632711038949265
No 21 is in faulty  (2,72,72) =  -0.0007396509740271954 (correct value)


-
- 72
- 1 2 
- - - -

N  = 3*24+1, Nd = 3*24-1

2187/246400


=#
