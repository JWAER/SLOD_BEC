function Coarse2Fine(Mesh)

P_h = spzeros(size(Mesh.p,2),size(Mesh.p,2));

I = ones(Int,10^8);
J = ones(Int,10^8);

idx_IJ = 1;
for i = 1:size(Mesh.t,2);

	Coarse_Nodes = Mesh.t[[1,4,10,20],i];
	Coarse_Points = Mesh.p[:,Coarse_Nodes];
	


        Jacb = hcat(Coarse_Points[:,2]-Coarse_Points[:,1], Coarse_Points[:,3]-Coarse_Points[:,1], Coarse_Points[:,4]-Coarse_Points[:,1]);


        DetJ = det(Jacb);

	x_0 = Coarse_Points[:,1];
	idx_f = Mesh.t[:,i];
	idx_C = Mesh.t[[1,4,10,20],i];

	for j = 1 #fine simplices
		for k = 1:20
			for l = 1:4
	 			I[idx_IJ] = idx_C[l]; J[idx_IJ] = idx_f[k]; idx_IJ +=1;#P_h[idx_C[l],idx_f[k]]= abs(Vals[l]);
			end
			

		end
	end


	if(idx_IJ>0.9*10^8); P_h+= sparse(I[1:idx_IJ-1],J[1:idx_IJ-1],2^-100,size(Mesh.p,2),size(Mesh.p,2)); fill!(I,1);fill!(J,1); idx_IJ = 1; end
	
end


	P_h += sparse(I[1:idx_IJ-1],J[1:idx_IJ-1],2^-100,size(Mesh.p,2),size(Mesh.p,2));
	
	
for i = 1:size(Mesh.t,2);

	Coarse_Nodes = Mesh.t[[1,4,10,20],i];
	Coarse_Points = Mesh.p[:,Coarse_Nodes];
	


        Jacb = hcat(Coarse_Points[:,2]-Coarse_Points[:,1], Coarse_Points[:,3]-Coarse_Points[:,1], Coarse_Points[:,4]-Coarse_Points[:,1]);


        DetJ = det(Jacb);

	x_0 = Coarse_Points[:,1];
	idx_f = Mesh.t[:,i];
	idx_C = Mesh.t[[1,4,10,20],i];

	#for j = 1 #fine simplices
		for k = 1:20
			point = Mesh.p[:,Mesh.t[k,i]];
	
			X = Jacb\(point-x_0)
			Vals = [1-sum(X),X[1],X[2],X[3]]
			for l = 1:4
	 			P_h[idx_C[l],idx_f[k]]= abs(Vals[l]);
			end
			

		end
	#end

	
end


	
	

return P_h;
end







function Assemble_P(Mesh)
	Dim = size(Mesh.p,2);

	I_vec = zeros(Int,10^6);
	J_vec = zeros(Int,10^6);
	idx = 1;
	
	P = spzeros(Dim,Dim);	

	
	#-----------------------------------Allocate------------------------------------------------
	for s = 1:size(Mesh.t_H,2);
		dofs_C = Mesh.t_H[:,s];
		
		fine_simplices = Mesh.FineSimplicesInCoarse[s];
		dofs_f = unique(Mesh.t[:,fine_simplices])
		len = length(dofs_f);
		
		for i = 1:3
			I_vec[idx:idx+len-1] .= dofs_C[i]; J_vec[idx:idx+len-1] = dofs_f; idx += len;
		end
		
		
		if(idx>0.8*10^6); P+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],eps(10^-10),Dim,Dim); idx = 1; end
		
	end
	
	P+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],eps(10^-10),Dim,Dim);
	
	
	
	
	#---------------------------------COMPUTE -----------------------------------"

	for s = 1:size(Mesh.t_H,2);
	
		P1_nodes = Mesh.t_H[:,s];
		
		vertices = Mesh.p[:,P1_nodes]
	        J = vertices[:,2:4].-vertices[:,1];
		DetJ = det(J);
		
		fine_simplices = Mesh.FineSimplicesInCoarse[s];
		nodes_f = unique(Mesh.t[:,fine_simplices])
		len = length(nodes_f);
	
		for node in nodes_f
			p_n = Mesh.p[:,node];
			X = J\(p_n-vertices[:,1])
			P1_VALS = [1-sum(X),X[1],X[2],X[3]];;
			P1_VALS = abs.(P1_VALS.*(P1_VALS.>0));
			for i = 1:4
				P[P1_nodes[i],node] =P1_VALS[i];
			end
		end
		
		
	end
	
	return P
end

