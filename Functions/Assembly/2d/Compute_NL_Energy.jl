function NL_Energy(U,Quad,Mesh);

	N_gp = Int(Quad["N_gp"][1]);
	val = 0;
	

	for s = 1:size(Mesh.t_h,2);
		loc2glob = Mesh.loc2glob[:,s];
		Cs  = U[loc2glob]; #U on full space
		
		vertices = Mesh.p[:,Mesh.t_h[:,s]]
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1]);
		DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp];
		
		rho = 0; for ii = 1:10; rho+=vₕ[ii](X,Y)*Cs[ii]; end; rho = abs(rho)^2;
		val += rho^2*w*DetJ; 

		end
			
			
	end

	return val;


end
