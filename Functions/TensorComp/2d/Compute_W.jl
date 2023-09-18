# W[i][j,k] beräknas för j,k>= i vilket betyder att W[i][j,k] för j,k<i, antag k<=j<i, sätt w[i][j,k] = w[k][j,i]

#For each coarse simplex find quadrature points, then find corresponding fine simplex and do as usual!
#----------------------------------------------------------------------------------------------#

function  Compute(mesh,W,ϕ,dofs_f)

	#----------------Quadrature rule to EXACTLY integrate a 9th degree polynomial-----------------

X_gp = [ 
   0.020634961602525
   0.036838412054736
   0.036838412054736
   0.044729513394453
   0.044729513394453
   0.125820817014127
   0.188203535619033
   0.188203535619033
   0.221962989160766
   0.221962989160766
   0.333333333333333
   0.437089591492937
   0.437089591492937
   0.489682519198738
   0.489682519198738
   0.623592928761935
   0.741198598784498
   0.741198598784498
   0.910540973211095];
	   
	Y_gp = [   0.489682519198738
   0.221962989160766
   0.741198598784498
   0.044729513394453
   0.910540973211095
   0.437089591492937
   0.188203535619033
   0.623592928761935
   0.036838412054736
   0.741198598784498
   0.333333333333333
   0.125820817014127
   0.437089591492937
   0.020634961602525
   0.489682519198738
   0.188203535619033
   0.036838412054736
   0.221962989160766
   0.044729513394453];

	w_gp =  [   0.015667350113570
   0.021641769688645
   0.021641769688645
   0.012788837829349
   0.012788837829349
   0.038913770502387
   0.039823869463605
   0.039823869463605
   0.021641769688645
   0.021641769688645
   0.048567898141400
   0.038913770502387
   0.038913770502387
   0.015667350113570
   0.015667350113570
   0.039823869463605
   0.021641769688645
   0.021641769688645
   0.012788837829349];

	N_gp = length(w_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	w_gp = reshape(w_gp,N_gp,1);
	#-----------------------------------------------------------------------------------------   

	max_number = 100 #max number of basis functions with support on simplex REPLACE W. FORMULA
	quad_vals = rand(max_number,N_gp); #Contains function values in quadrature point



map_f = zeros(Int,size(mesh.p,2));
for i = 1:length(dofs_f); map_f[dofs_f[i]] = i; end


	ta = time();

	
#tid = time();

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)

		Simplex = mesh.t_h[:,n_t];
		Vertices =  mesh.p[:,Simplex];
		J = hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1]);
		Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp]; 
		#loc2glob = setdiff(map_f[mesh.loc2glob[:,n_t]],0);
		loc2glob = map_f[mesh.loc2glob[:,n_t]];
		
	OD_S = ϕ[:,setdiff(loc2glob,0)];
	ODval = Matrix(OD_S);
	ODs = similar(OD_S.rowval); copyto!(ODs,OD_S.rowval);
	unique!(sort!(ODs));
	len = length(ODs);
	
	#add to global tensor using local tensor which is constant!
	detJ = det(J)
	fill!(quad_vals,0.0)
	for gp in 1:N_gp#Quad2FineSimplex
	
	
	for ii = 1:10;if(loc2glob[ii]!=0); quad_vals[1:len,gp] += ϕ[ODs,loc2glob[ii]]*vₕ[ii](X_gp[gp],Y_gp[gp]); end; end
		
	end
		
	w_H = detJ*w_gp;
 	Quad_Rule(quad_vals,w_H,len,ODs,W)	
		
	end 


end


function Quad_Rule(quad_vals,w_H,len,a_n,W)

	for i = 1:len
		i_glob = a_n[i];

		sub_i = W.Iptr[i_glob]: W.Iptr[i_glob+1]-1;
		
if(!isempty(sub_i))

		PreAlloc_IdxJ = W.J[sub_i];

		f = quad_vals[i,:].*w_H;
			
		for j = i:len
			j_glob = a_n[j]
			sub_ij= (sub_i[1]-1).+searchsorted(PreAlloc_IdxJ, j_glob)

	
			g = f.*quad_vals[j,:];
			
		if(!isempty(sub_ij))	
	
			PreAlloc_IdxK = W.K[sub_ij]
			for k = j:len
				k_glob = a_n[k];
				sub_ijk=(sub_ij[1]-1).+searchsorted(PreAlloc_IdxK, k_glob)
				if(!isempty(sub_ijk)) W.Val[sub_ijk[1]]+=dot(quad_vals[k,:],g);  end
				#W.Val[sub_ijk].+=val
			end
		end
		end

	end
end

end




























#----------------------------------------------------------------------------------------#

function  Compute_Canonical(mesh,W,ϕ,dofs_f,ℓ, compute_at_bdry = false )

	#----------------Quadrature rule to EXACTLY integrate a 9th degree polynomial-----------------

X_gp = [ 
   0.020634961602525
   0.036838412054736
   0.036838412054736
   0.044729513394453
   0.044729513394453
   0.125820817014127
   0.188203535619033
   0.188203535619033
   0.221962989160766
   0.221962989160766
   0.333333333333333
   0.437089591492937
   0.437089591492937
   0.489682519198738
   0.489682519198738
   0.623592928761935
   0.741198598784498
   0.741198598784498
   0.910540973211095];
	   
	Y_gp = [   0.489682519198738
   0.221962989160766
   0.741198598784498
   0.044729513394453
   0.910540973211095
   0.437089591492937
   0.188203535619033
   0.623592928761935
   0.036838412054736
   0.741198598784498
   0.333333333333333
   0.125820817014127
   0.437089591492937
   0.020634961602525
   0.489682519198738
   0.188203535619033
   0.036838412054736
   0.221962989160766
   0.044729513394453];

	w_gp =  [   0.015667350113570
   0.021641769688645
   0.021641769688645
   0.012788837829349
   0.012788837829349
   0.038913770502387
   0.039823869463605
   0.039823869463605
   0.021641769688645
   0.021641769688645
   0.048567898141400
   0.038913770502387
   0.038913770502387
   0.015667350113570
   0.015667350113570
   0.039823869463605
   0.021641769688645
   0.021641769688645
   0.012788837829349];

	N_gp = length(w_gp)
	X_gp = reshape(X_gp,1,N_gp);
	Y_gp = reshape(Y_gp,1,N_gp);
	w_gp = reshape(w_gp,N_gp,1);
	#-----------------------------------------------------------------------------------------   

	max_number = 100 #max number of basis functions with support on simplex REPLACE W. FORMULA
	quad_vals = rand(max_number,N_gp); #Contains function values in quadrature point

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);
map_C = [searchsorted(mesh.CoarseNodes,dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that 


center_node_C = find_center(dofs_C);

map_f = zeros(Int,size(mesh.p,2));
for i = 1:length(dofs_f); map_f[dofs_f[i]] = i; end


	ta = time();

	
#tid = time();

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)

	Simplex = mesh.t_h[:,n_t];
			
	Vertices =  mesh.p[:,Simplex];
		
	dist2bdry = min( minimum(abs.(Vertices.-mesh.box)), minimum(abs.(Vertices.+mesh.box)))
	
	
		J = hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1]);
		Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp]; 
		#loc2glob = setdiff(map_f[mesh.loc2glob[:,n_t]],0);
		loc2glob = map_f[mesh.loc2glob[:,n_t]];
		
		OD_S = ϕ[:,setdiff(loc2glob,0)];
		ODval = Matrix(OD_S);
		ODs = similar(OD_S.rowval); copyto!(ODs,OD_S.rowval);
		unique!(sort!(ODs));
		len = length(ODs);
		
		#if( (center_node_C in ODs ) | (dist2bdry < (4ℓ+2)*mesh.H)*compute_at_bdry );	
		if( (center_node_C in ODs ) );	
	
	
		
		#add to global tensor using local tensor which is constant!
		detJ = det(J)
		fill!(quad_vals,0.0)
		for gp in 1:N_gp#Quad2FineSimplex
		
		
		for ii = 1:10;if(loc2glob[ii]!=0); quad_vals[1:len,gp] += ϕ[ODs,loc2glob[ii]]*vₕ[ii](X_gp[gp],Y_gp[gp]); end; end
			
		end
			
		w_H = detJ*w_gp;
	 	Quad_Rule(quad_vals,w_H,len,ODs,W)	
		end
	
	
	end
	idx = W.Iptr[center_node_C]:W.Iptr[center_node_C+1]-1;
	standard_values = W.Val[idx];
	
	len = length(idx);
	
	for i = 1:length(W.Iptr)-1  #copy standard values (and overwrite certain) 
	
		point = mesh.p[:,mesh.CoarseNodes[map_C[i]]];
		dist2bdry = min( minimum(abs.(point .-mesh.box)), minimum(abs.(point.+mesh.box)))
	
	#	if( ((W.Iptr[i+1]-W.Iptr[i]) == len) && (dist2bdry> ( (3ℓ+2)*mesh.H)) );
		if( ((W.Iptr[i+1]-W.Iptr[i]) == len) );
	
			W.Val[W.Iptr[i]:W.Iptr[i+1]-1].=standard_values;
		end	
	end
	

end






function  Compute_Wh(mesh,W)


	#----------------Quadrature rule to EXACTLY integrate a 9th degree polynomial-----------------

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

	active_nodes = Vector{Int64}(undef,2000); #4*300 is the underlying thought 
	W_loc = zeros(10,10,10);

dofs_f = setdiff(1:size(mesh.p,2),mesh.bdry);
map_f = zeros(Int,size(mesh.p,2));
for i = 1:length(map_f); if(!isempty(searchsorted(dofs_f,i)));map_f[i] = searchsorted(dofs_f,i)[1]; end; end

loc2glob_dirichlet = map_f[mesh.loc2glob];

	ta = time();
	
#tid = time();

	for n_t = 1:size(mesh.t_h,2)#size(t_C,2)
	
		Simplex = mesh.t_h[:,n_t];
		Vertices =  mesh.p[:,Simplex];
		J = hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1]);
		Quad_Points = Vertices[:,1].+J*[X_gp;Y_gp]; 
		loc2glob =loc2glob_dirichlet[:,n_t]
	
	
	if(n_t==1)	
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
	
	#add to global tensor using local tensor which is constant!
	detJ = det(J)

	Quad_Rule(W_loc,detJ,loc2glob,W) 	
	
		
	end 


end




function Quad_Rule(W_loc,detJ,loc2glob,W)
#sort loc2glob? or use if-statements

	for i = 1:10
		i_glob = loc2glob[i];
	if(i_glob!=0)
		sub_i = W.Iptr[i_glob]: W.Iptr[i_glob+1]-1;

		PreAlloc_IdxJ = W.J[sub_i];

		for j = 1:10
			j_glob = loc2glob[j]
		if(i_glob<=j_glob);	
			sub_ij= (sub_i[1]-1).+searchsorted(PreAlloc_IdxJ, j_glob)

		if(!isempty(sub_ij))	
	
			PreAlloc_IdxK = W.K[sub_ij]
			for k = 1:10
			k_glob = loc2glob[k];
			if(j_glob<=k_glob)
				sub_ijk=(sub_ij[1]-1).+searchsorted(PreAlloc_IdxK, k_glob)
				if(!isempty(sub_ijk)) W.Val[sub_ijk[1]]+=W_loc[i,j,k]*detJ  end
			end
			end
		end
		end
		end

	end
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






function remove_zeros(W_new,W,SpaceDim_C)
	count = 1;
	for i = 1:SpaceDim_C
		for idx = W.Iptr[i]:W.Iptr[i+1]-1
			if(W_val[idx] == 0)
				W_new.Iptr[i+1:end] .-=1;
			else
				W_new.J[count] = W.J[idx];
				W_new.K[count] = W.K[idx];
				W_new.val[count] = W_val[idx];	
				count += 1;
			end
		end
	end

end




function find_center(dofs_C)

	N = Int(sqrt(length(dofs_C)));
	
	if(mod(N,2)!=0); return Int(round(N*N/2)); end
	
	return Int( N*N/2-N/2);
end























