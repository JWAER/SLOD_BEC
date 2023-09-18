function Compute_RHS(V,t,Quad,Mesh,dofs_f);

	rhs = zeros(size(Mesh.p,2));

	N_gp = Int(Quad["N_gp"][1]);
	val = 0;
	
	
	pre_alloc = zeros(10,N_gp);
	for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		for ii = 1:10; pre_alloc[ii,gp]=vₕ[ii](X,Y) end; 
	end
		
	
	
	

	for s = 1:size(Mesh.t_h,2);
		loc2glob = Mesh.loc2glob[:,s];
		
		vertices = Mesh.p[:,Mesh.t_h[:,s]]
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1]);
		DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp];
		x = vertices[:,1]+J*[X;Y];
		
		val = V(x,t)*w*DetJ;;
		
		for ii = 1:10; rhs[loc2glob[ii]]+=pre_alloc[ii,gp]*val;end#ϕ[ii](X,Y)*Vval; end; 
		#for ii = 1:10; rhs[loc2glob[ii]]+=ϕ[ii](X,Y)*val; end; 
	
		end
			
			
	end

	return rhs[dofs_f];


end
