# W[i][j,k] beräknas för j,k>= i vilket betyder att W[i][j,k] för j,k<i, antag k<=j<i, sätt w[i][j,k] = w[k][j,i]

#For each coarse simplex find quadrature points, then find corresponding fine simplex and do as usual!
#----------------------------------------------------------------------------------------------#

function  Compute_W(mesh, ω,ϕ,Quad)
	#----------------Quadrature rule to EXACTLY integrate a 9th degree polynomial-----------------

	X_gp = Quad["X"];
	Y_gp = Quad["Y"];
	Z_gp = Quad["Z"];
	 ω_gp = Quad["ω"];

	N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	Z_gp = reshape(Z_gp,1,N_gp);
	ω_gp = reshape(ω_gp,N_gp,1);
#=
	X_gp = load("XGP.jld")["XGP"];
	Y_gp = load("YGP.jld")["YGP"];
	Z_gp = load("ZGP.jld")["ZGP"];
	ω_gp = load("WGP.jld")["WGP"];

	N_gp = length(ω_gp)
		X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	Z_gp = reshape(Y_gp,1,N_gp);
	ω_gp = reshape(ω_gp,N_gp,1)/27;
	
	println(size(ω_gp))
=#
	#-----------------------------------------------------------------------------------------   

	max_number = 100 #max number of basis functions with support on simplex REPLACE W. FORMULA
	quad_vals = rand(max_number,N_gp); #Contains function values in quadrature point


map_f = zeros(Int,size(mesh.p,2));
for i = 1:length(mesh.dofs); map_f[mesh.dofs[i]] = i; end


	ta = time();

	
#tid = time();

	for n_t = 1:size(mesh.t,2)#size(t_C,2)
	

		Simplex = mesh.t[[1,4,10,20],n_t];
		Vertices =  mesh.p[:,Simplex];
		J = Vertices[:,2:end].-Vertices[:,1]
	
		
	#	Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp;Z_gp]; 
		#loc2glob = setdiff(map_f[mesh.loc2glob[:,n_t]],0);
		loc2glob = map_f[mesh.t[:,n_t]];
		
	LOD_S = ϕ[:,setdiff(loc2glob,0)];
	LODval = Matrix(LOD_S);
	LODs = similar(LOD_S.rowval); copyto!(LODs,LOD_S.rowval);
	unique!(sort!(LODs)); #These have support on simplex
	len = length(LODs);
	 
	#add to global tensor using local tensor which is constant!
	detJ = det(J)

	fill!(quad_vals,0.0)
	for gp in 1:N_gp#Quad2FineSimplex
	
	
	for ii = 1:20;if(loc2glob[ii]!=0); quad_vals[1:len,gp] += ϕ[LODs,loc2glob[ii]]*vₕ[ii]([X_gp[gp],Y_gp[gp],Z_gp[gp]]); end; end


		
	end
	w_H = detJ*ω_gp;

 	Quad_Rule(quad_vals,w_H,len,LODs, ω)	
		

	end 


end


function Quad_Rule(quad_vals,w_H,len,a_n, ω)

	for i = 1:len
		i_glob = a_n[i];

		sub_i =  ω.Iptr[i_glob]:  ω.Iptr[i_glob+1]-1;
		
if(!isempty(sub_i))

		PreAlloc_IdxJ =  ω.J[sub_i];

		f = quad_vals[i,:].*w_H;
			
		for j = i:len
			j_glob = a_n[j]
			sub_ij= (sub_i[1]-1).+searchsorted(PreAlloc_IdxJ, j_glob)

	
			g = f.*quad_vals[j,:];
			
		if(!isempty(sub_ij))	
	
			PreAlloc_IdxK = ω.K[sub_ij]
			for k = j:len
				k_glob = a_n[k];
				sub_ijk=(sub_ij[1]-1).+searchsorted(PreAlloc_IdxK, k_glob)
				if(!isempty(sub_ijk)) ω.Val[sub_ijk[1]]+=dot(quad_vals[k,:],g);  end
				
			end
		end
		end

	end
end

end









# W[i][j,k] beräknas för j,k>= i vilket betyder att W[i][j,k] för j,k<i, antag k<=j<i, sätt w[i][j,k] = w[k][j,i]

#For each coarse simplex find quadrature points, then find corresponding fine simplex and do as usual!
#----------------------------------------------------------------------------------------------#

function  Compute_Canonical(mesh, ω,ϕ,Quad,ℓ,compute_at_bdry = false)


	X_gp = Quad["X"];
	Y_gp = Quad["Y"];
	Z_gp = Quad["Z"];
	 ω_gp = Quad["ω"];

	N_gp = length(ω_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	Z_gp = reshape(Z_gp,1,N_gp);
	ω_gp = reshape(ω_gp,N_gp,1);

	#-----------------------------------------------------------------------------------------   

	max_number = 200 #max number of basis functions with support on simplex REPLACE W. FORMULA
	quad_vals = rand(max_number,N_gp); #Contains function values in quadrature point


map_f = zeros(Int,size(mesh.p,2));
for i = 1:length(mesh.dofs); map_f[mesh.dofs[i]] = i; end

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);
map_C = [searchsorted(mesh.CoarseNodes,dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C]

center_node_C = Int(round(length(map_C)/2+0.01));
center_node =mesh.CoarseNodes[ map_C[center_node_C]];#Compute this one first then translate, whenever possible, i.e.,. when lengths are equal


	ta = time();

	
#tid = time();

	for n_t = 1:size(mesh.t,2)#size(t_C,2)
	

		Simplex = mesh.t[[1,4,10,20],n_t];
		Vertices =  mesh.p[:,Simplex];
		J = Vertices[:,2:end].-Vertices[:,1]
		
		dist2bdry = min(minimum(abs.(Vertices.-mesh.box_size)), minimum(abs.(Vertices.+mesh.box_size)))
		
	#	Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp;Z_gp]; 
		#loc2glob = setdiff(map_f[mesh.loc2glob[:,n_t]],0);
		loc2glob = map_f[mesh.t[:,n_t]];
		
	LOD_S = ϕ[:,setdiff(loc2glob,0)];
	LODval = Matrix(LOD_S);
	LODs = similar(LOD_S.rowval); copyto!(LODs,LOD_S.rowval);
	unique!(sort!(LODs)); #These have support on simplex
	len = length(LODs);
	
	if( (center_node_C in LODs ) | (dist2bdry < (4ℓ+2)*mesh.H)*compute_at_bdry );	
	
		 
		#add to global tensor using local tensor which is constant!
		detJ = det(J)

		fill!(quad_vals,0.0)
		for gp in 1:N_gp#Quad2FineSimplex
		
		
			for ii = 1:20;if(loc2glob[ii]!=0); quad_vals[1:len,gp] += ϕ[LODs,loc2glob[ii]]*vₕ[ii]([X_gp[gp],Y_gp[gp],Z_gp[gp]]); end; end


				
			end
			w_H = detJ*ω_gp;

		 	Quad_Rule(quad_vals,w_H,len,LODs, ω)	
		end

	end 
	
	idx =  ω.Iptr[center_node_C]: ω.Iptr[center_node_C+1]-1;
	standard_values = ω.Val[idx];
	
	len = length(idx);
	
	for i = 1:length( ω.Iptr)-1  #copy standard values (and overwrite certain) 
		
		point = mesh.p[:,mesh.CoarseNodes[map_C[i]]];
		dist2bdry = min( minimum(abs.(point .-mesh.box_size)), minimum(abs.(point.+mesh.box_size)))
		
		
		if( (( ω.Iptr[i+1]- ω.Iptr[i]) == len) );
			println(i);
			ω.Val[ ω.Iptr[i]: ω.Iptr[i+1]-1].=standard_values;
		end
#=
	
		if( (( ω.Iptr[i+1]- ω.Iptr[i]) == len) && (dist2bdry>= ( (3ℓ+2)*mesh.H)-eps(10.0) ));
			println(i);
			ω.Val[ ω.Iptr[i]: ω.Iptr[i+1]-1].=standard_values;
		
		
		end
=#	
	end
	


end

function Compute_tilde(W)
prealloc = length(W.Val);
I_vec = zeros(Int,prealloc);
J_vec = zeros(Int,prealloc);
K_vec = zeros(Int,prealloc);
V_vec = zeros(prealloc);
it=1;

N = length(W.Iptr)-1;



for i = 1:length(W.Iptr)-1
	to_do = W.Iptr[i]:W.Iptr[i+1]-1;

	for idx = to_do;
		j = W.J[idx];
		k = W.K[idx];
		if(i!=k && i!=j)
		I_vec[it] = j; J_vec[it] = k; K_vec[it] = i; V_vec[it]= W.Val[idx];it+=1;
		end		
	end		
end
	
I_vec = I_vec[1:it-1]; J_vec = J_vec[1:it-1]; K_vec = K_vec[1:it-1]; V_vec = V_vec[1:it-1];

p = sortperm(I_vec);

I_vec = I_vec[p]; J_vec = J_vec[p]; K_vec = K_vec[p]; V_vec = V_vec[p];

Iptr = zeros(Int,N+1);
Iptr[1] = 1;



for i = 1:N
	sub_i = searchsorted(I_vec,i);
	
	Iptr[i+1] = Iptr[i]+length(sub_i); 

	p_ = sortperm(J_vec[sub_i]);
	J_vec[sub_i] = J_vec[sub_i][p_];
	K_vec[sub_i] = K_vec[sub_i][p_];
	V_vec[sub_i] = V_vec[sub_i][p_];

end


return Tensor(Iptr,J_vec,K_vec,V_vec);

end




