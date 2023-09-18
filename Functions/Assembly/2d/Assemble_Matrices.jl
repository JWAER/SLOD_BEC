function MatrixSparsity(mesh)
	N = size(mesh.p,2);
	M = spzeros(N,N);
	
	Idx_I = ones(Int,10^7);
   	Idx_J = ones(Int,10^7);
   	idx_IJ = 1;
   	
   	loc2glob = zeros(Int,10);
   
for t = 1:size(mesh.t_h,2);


        loc2glob .= mesh.loc2glob[:,t];

	    for i = 1:10
	       i_glob = loc2glob[i];
	        for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
	            Idx_I[idx_IJ] = i_glob; Idx_J[idx_IJ] = j_glob;  idx_IJ += 1;	     
	        end end
	    end
	



  
     #---Local matrix calculated------

     #Added   
	if(idx_IJ>0.8*10^7) M+= sparse(Idx_I[1:idx_IJ-1],Idx_J[1:idx_IJ-1],1.0,N,N); fill!(Idx_I,1);fill!(Idx_J,1);idx_IJ = 1; end
	

end

	M+= sparse(Idx_I[1:idx_IJ-1],Idx_J[1:idx_IJ-1],1.0,N,N); 
	
	M+=M'; 
	
	fill!(M.nzval,0.)

return M


end





function vᵢvⱼ_OLD(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);
	
	M = spzeros(Dim,Dim);	copyto!(M,sparsity);
	M_loc = zeros(10,10);
	
	vertices = zeros(2,3);
	J = zeros(2,2);
	loc2glob = zeros(Int,10);	
		

	for s = 1:size(Mesh.t_h,2);
		fill!(M_loc,0.);
		loc2glob .= Mesh.loc2glob[:,s];
	
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J .= vertices[:,2:3].-vertices[:,1];
		DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					M_loc[j,i]+= loc_i*vₕ[j](X,Y)*w*DetJ; 
				end
			end 
		end
		M_loc+=M_loc'; for i = 1:10; M_loc[i,i]/=2; end
		
	
		for i = 1:10
		       i_glob = loc2glob[i];
		       idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		       row = M.rowval[idx];
		
		        for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		        m = searchsorted(row,j_glob)[1];
			M.nzval[idx[m]]+=M_loc[j,i];
				
		       end end
		end
	
		
	end
	
	M+=M'; for d = 1:size(M,1); M[d,d]/=2; end
	

	return M
end



function vᵢvⱼ(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);
	
	M = spzeros(Dim,Dim);	copyto!(M,sparsity);
	M_loc = zeros(10,10);
	
	vertices = zeros(2,3);
	J = zeros(2,2);
	loc2glob = zeros(Int,10);	
		
for modular = 1:2;

	for s = modular:2:size(Mesh.t_h,2);
	loc2glob .= Mesh.loc2glob[:,s];
	
	
	if(s==modular)
		fill!(M_loc,0.);
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J .= vertices[:,2:3].-vertices[:,1];
		DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					M_loc[j,i]+= loc_i*vₕ[j](X,Y)*w*DetJ; 
				end
			end 
		end
		M_loc+=M_loc'; for i = 1:10; M_loc[i,i]/=2; end
	end	
	
		for i = 1:10
		       i_glob = loc2glob[i];
		       idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		       row = M.rowval[idx];
		
		        for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		        m = searchsorted(row,j_glob)[1];
			M.nzval[idx[m]]+=M_loc[j,i];
				
		       end end
		end
	
		
	end
	
end
	M+=M'; for d = 1:size(M,1); M[d,d]/=2; end
	

	return M
end


#----------------------------------------------------------------------------------------------------------------------------------
function Vvᵢvⱼ(V,Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);
	
	
	M = spzeros(Dim,Dim);	copyto!(M,sparsity)
	M_loc = zeros(10,10);
	
	loc2glob = zeros(Int,10);
	vertices = zeros(2,3);
	J = zeros(2,2);
	

	for s = 1:size(Mesh.t_h,2);
		fill!(M_loc,0.0);
		loc2glob .= Mesh.loc2glob[:,s];
	
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J = vertices[:,2:3].-vertices[:,1];
	        DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
		
		x = vertices[:,1]+J*[X;Y];
		
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					M_loc[j,i] += loc_i*vₕ[j](X,Y)*V(x)*w*DetJ; 
				end
			end 
		end
		
		M_loc +=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
		for i = 1:10
			i_glob = loc2glob[i];
			idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
			row = M.rowval[idx];
		
			for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
			m = searchsorted(row,j_glob)[1];
			M.nzval[idx[m]]+=M_loc[j,i];
				
		       end end
		end
				
	end
	
	M +=M'; for d = 1:size(M,1); M[d,d]/=2; end
	
	return M
end

#-------------------------------------------------------------------------------------------------------------


function ∇vᵢ∇vⱼ_OLD(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);

	A = spzeros(Dim,Dim);	copyto!(A,sparsity);
	A_loc = zeros(10,10);
	
	loc2glob = zeros(Int,10);
	vertices = zeros(2,3);
	J = zeros(2,2);
	

	for s = 1:size(Mesh.t_h,2);
		fill!(A_loc,0);
		
		loc2glob .= Mesh.loc2glob[:,s];
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J = vertices[:,2:3].-vertices[:,1];
	        DetJ = det(J);
		JTinv = inv(J)';	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
			for i=1:10
				grad_i = JTinv*∇vₕ[i](X,Y);
				for j = i:10
				grad_j = JTinv*∇vₕ[j](X,Y);
					A_loc[j,i]+= dot(grad_i,grad_j)*w*DetJ; 
				end
			end 
		end
		
		A_loc+=A_loc'; for d = 1:10 A_loc[d,d]/=2; end
		
		for i = 1:10
			i_glob = loc2glob[i];
			idx = A.colptr[i_glob]:A.colptr[i_glob+1]-1;
			row = A.rowval[idx];
		
			for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
			m = searchsorted(row,j_glob)[1];
			A.nzval[idx[m]]+=A_loc[j,i];
				
		       end end
		end
					
		
	end
	
	A += A'; for d = 1:size(A,1); A[d,d]/=2; end
	
	return A
end





function ∇vᵢ∇vⱼ(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);

	A = spzeros(Dim,Dim);	copyto!(A,sparsity);
	A_loc = zeros(10,10);
	
	loc2glob = zeros(Int,10);
	vertices = zeros(2,3);
	J = zeros(2,2);


for modular  = 1:2	

	for s = modular:2:size(Mesh.t_h,2);
		
		
		loc2glob .= Mesh.loc2glob[:,s];
		
	if(s==modular) #compute reference A_loc
		fill!(A_loc,0);
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J = vertices[:,2:3].-vertices[:,1];
	        DetJ = det(J);
		JTinv = inv(J)';	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
			for i=1:10
				grad_i = JTinv*∇vₕ[i](X,Y);
				for j = i:10
				grad_j = JTinv*∇vₕ[j](X,Y);
					A_loc[j,i]+= dot(grad_i,grad_j)*w*DetJ; 
				end
			end 
		end
		
		A_loc+=A_loc'; for d = 1:10 A_loc[d,d]/=2; end
	end	
		for i = 1:10
			i_glob = loc2glob[i];
			idx = A.colptr[i_glob]:A.colptr[i_glob+1]-1;
			row = A.rowval[idx];
		
			for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
			m = searchsorted(row,j_glob)[1];
			A.nzval[idx[m]]+=A_loc[j,i];
				
		       end end
		end
					
		
	end
	
end
	
	A += A'; for d = 1:size(A,1); A[d,d]/=2; end
	
	return A
end














function P(Mesh)
	Dim = size(Mesh.p,2);

	I_vec = zeros(Int,10^6);
	J_vec = zeros(Int,10^6);
	idx = 1;
	
	P = spzeros(Dim,Dim);	

	
	#-----------------------------------Allocate------------------------------------------------
	for s = 1:size(Mesh.t_H,2);
		dofs_C = Mesh.t_H[:,s];
		
		fine_simplices = Mesh.FineSimplicesInCoarse[s];
		dofs_f = unique(Mesh.loc2glob[:,fine_simplices])
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
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1]);
		DetJ = det(J);
		
		fine_simplices = Mesh.FineSimplicesInCoarse[s];
		nodes_f = unique(Mesh.loc2glob[:,fine_simplices])
		len = length(nodes_f);
	
		for node in nodes_f
			p_n = Mesh.p[:,node];
			X = J\(p_n-vertices[:,1])
			P1_VALS = [1-sum(X),X[1],X[2]];
			P1_VALS = abs.(P1_VALS.*(P1_VALS.>0));
			for i = 1:3
				P[P1_nodes[i],node] =P1_VALS[i];
			end
		end
		
		
	end
	
	return P
end






#--------------------------------------------------------------------#
#COMPUTE <v_i,d_n(v_j)>

function M_bc(Mesh,bd_f,bd);
bd_tot = zeros(Int,3*length(bd));
idx_bd = 1;

Nbc = length(bd);
X_gp = [1 1+sqrt(3/5) 1-sqrt(3/5)]/2;
w = [4/9 5/18 5/18];

I_vec = zeros(Int,10^5);
J_vec = zeros(Int,10^5);
V_vec = zeros(10^5);
it = 1;
	for i = 1:Nbc

		possible_simplices= intersect(Mesh.FineNode2FineSimplices[:,bd[i]],Mesh.FineNode2FineSimplices[:,bd[mod(i,Nbc)+1]]);
		setdiff!(possible_simplices,0);
	
		s1 = possible_simplices[1];
		n1 = bd[i]; n2= bd[mod(i,Nbc)+1]; n3 = setdiff(Mesh.t_h[:,s1],[n1,n2])[1];
	
		test=det(Mesh.p[:,[n2,n3]].-Mesh.p[:,n1])

		if(test>0);s = s1; else s = possible_simplices[2];end

		idx = findfirst(Mesh.t_h[:,s].==n1)[1]; #LOCATES IDX OF FIRST BDRY NODE, Right hand rule implies next
		
		loc2glob = Mesh.loc2glob[:,s];	
		vertices = Mesh.p[:,Mesh.t_h[:,s]]
	
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1]);
		JTinv = inv(J)';
			
		p1 = vertices[:,idx]; 		#first bdnode
		p2 = vertices[:,mod(idx,3)+1] #second bdnode
		
		P1 = [0 1 0 ; 0 0 1][:,idx];
		P2 = [0 1 0 ; 0 0 1][:,mod(idx,3)+1];
		
		H = norm(p2-p1);

		normal = [0 1; -1 0]*(p2-p1); normalize!(normal);

		if(idx==1); support_on_bdry = [1,2,3,4];
		elseif(idx==2); support_on_bdry=[4,7,9,10];
		else; support_on_bdry=[1,5,8,10]; end
		
		for gp = 1:3
		x = p1+(p2-p1)*X_gp[gp]
		X = P1+(P2-P1)*X_gp[gp]
			for j = 1:10; #LOOP OVER ALL BASIS FUNCTIONS compute their normal deriv.
			gradn_j = dot(JTinv*∇vₕ[j](X[1],X[2]),normal)
			for ii = support_on_bdry#1:10; #Only 4 basis functions have support on bdry
				I_vec[it] = loc2glob[ii]; J_vec[it] = loc2glob[j]; V_vec[it]=gradn_j*vₕ[ii](X[1],X[2])*w[gp]*H; it+=1
			end
			end
		end
		
	end

dim = size(Mesh.p,2);


return sparse(I_vec[1:it-1],J_vec[1:it-1],V_vec[1:it-1],dim,dim)
end


function Vvᵢvⱼ(V,Mesh,Quad,Nh,sparsity)

ngp = length(Quad["X"]);
X_gp = zeros(Nh^2*ngp);
Y_gp = zeros(Nh^2*ngp);
w_gp = zeros(Nh^2*ngp);

N_gp = Nh^2*ngp;

#Refine reference simplex then add gps! x(X) = x0+(x1-x0)*X+(x2-x0)*Y
h_loc = 1/Nh;
it = 1;
for j = 0:Nh-1
	Y = j/Nh
	for i = 0:Nh-1-j
		X = i/Nh
		
		#do lower triangle 
		#x0_h = F(X,Y)
		#x1_h = F(X+h_loc,Y);
		#x2_h = F(X,Y+h_loc);
		X_gp[it:it+ngp-1] = Quad["X"]/Nh.+X;
		Y_gp[it:it+ngp-1] = Quad["Y"]/Nh.+Y;
		w_gp[it:it+ngp-1] = Quad["ω"]/Nh^2;
		it+=ngp	
			
		if(i<Nh-1-j) #DO UPPER
		X_gp[it:it+ngp-1] = X+h_loc.-Quad["X"]/Nh;
		Y_gp[it:it+ngp-1] = Y+h_loc.-Quad["Y"]/Nh;	
		w_gp[it:it+ngp-1] = Quad["ω"]/Nh^2;	
		it+=ngp	
		end
	end
end


	Dim = size(Mesh.p,2);
	
	I_vec = zeros(Int,10^6);
	J_vec = zeros(Int,10^6);
	Val_vec = zeros(10^6);
	idx = 1;
	
	M = spzeros(Dim,Dim);	

	

	for s = 1:size(Mesh.t_h,2);
		loc2glob = Mesh.loc2glob[:,s];
	
		vertices = Mesh.p[:,Mesh.t_h[:,s]]
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1]);
		DetJ = det(J);
	
		
		for gp = 1:N_gp
		X = X_gp[gp]
		Y = Y_gp[gp]
		w = w_gp[gp]
		x = vertices[:,1]+J*[X;Y];
		
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					I_vec[idx] = loc2glob[i]; J_vec[idx] = loc2glob[j]; Val_vec[idx] = loc_i*vₕ[j](X,Y)*V(x)*w*DetJ; idx+=1;
				end
			end 
		end
		
		
		if(idx>0.8*10^6); M+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],Val_vec[1:idx-1],Dim,Dim); idx = 1; end
		
	end
	
	M+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],Val_vec[1:idx-1],Dim,Dim);
	M+=M'; for i = 1:size(M,1); M[i,i]/=2.; end
	return M
end
















function M_Ω(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);

	M_Ω = 0im*spzeros(Dim,Dim);copyto!(M_Ω,sparsity);
	M_loc = zeros(10,10);
		
	loc2glob = zeros(Int,10);
	vertices = zeros(2,3);
	J = zeros(2,2);

	for s = 1:size(Mesh.t_h,2);
		fill!(M_loc,0.);
	
		loc2glob .= Mesh.loc2glob[:,s];
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J = vertices[:,2:3].-vertices[:,1];
		DetJ = det(J);
		JTinv = inv(J)';	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		w = Quad["ω"][gp]
		yx = reverse(vertices[:,1]+J*[X;Y]);
		yx[2]*=-1;
		
			for i=1:10
				grad_i = JTinv*∇vₕ[i](X,Y);
				for j = 1:10
				#grad_j = JTinv*∇vₕ[j](X,Y);
					M_loc[i,j] += -dot(yx,grad_i)*vₕ[j](X,Y)*w*DetJ;  #provided dirichlet bcs
				end
			end 
		end
		
		
		for i = 1:10
			i_glob = loc2glob[i];
			idx = M_Ω.colptr[i_glob]:M_Ω.colptr[i_glob+1]-1;
			row = M_Ω.rowval[idx];
		
			for j = 1:10;   j_glob = loc2glob[j]; 
			m = searchsorted(row,j_glob)[1];
			M_Ω.nzval[idx[m]]+=M_loc[i,j];
				
		       end
		end		
		

	end
	
	return M_Ω
end












#-------------Adapted to piecewise constant functions---------------------------------

function Vvᵢvⱼ_example2(V,Mesh,Quad,Nh,sparsity)
#Remove score parameter to make code work for any piecewise constant function.
#Currently tailored so that mod(x,3) == 0 lines are not refined (due to the mesh aligning with them)

function add_to_sparse(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)

	M_loc = zeros(10,10);

	for gp = 1:N_gp
		X = X_gp[gp]
		Y = Y_gp[gp]
		w = w_gp[gp]
		x = vertices[:,1]+J*[X;Y];
		
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					M_loc[j,i]+= loc_i*vₕ[j](X,Y)*V(x)*w*DetJ;
				end
			end 
	end
	
	M_loc+=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
	for i = 1:10
		i_glob = loc2glob[i];
		idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		row = M.rowval[idx];
	
		for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
		
	
end	




function add_to_sparse_REF(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)

	Jinv = 1/det(J)*[ J[2,2] -J[1,2]; -J[2,1] J[1,1]]
	M_loc = zeros(10,10);

	for gp = 1:N_gp
		x = [X_gp[gp]; Y_gp[gp]];
		w = w_gp[gp]
		
		#x = vertices[:,1]+J*[X;Y];
		X = Jinv*(x-vertices[:,1])
			for i=1:10
				loc_i=vₕ[i](X[1],X[2])
				for j = i:10
					M_loc[j,i]+=loc_i*vₕ[j](X[1],X[2])*V(x)*w; 
				end
			end 
	end
	
	M_loc +=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
	for i = 1:10
		i_glob = loc2glob[i];
		idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		row = M.rowval[idx];
	
		for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
		
	
end	



function refine(vertices);
points = [vertices (vertices[:,1]+vertices[:,2])/2 (vertices[:,2]+vertices[:,3])/2 (vertices[:,3]+vertices[:,1])/2];
t = [1,4,6,4,2,5,6,5,3,4,5,6];
return points[:,t];

end



X_gp = Quad["X"];
Y_gp = Quad["Y"];
w_gp = Quad["ω"];

N_gp = length(X_gp);


Xr_gp = zeros(10^6);
Yr_gp = zeros(10^6);
wr_gp = zeros(10^6);

Dim = size(Mesh.p,2);

M = spzeros(Dim,Dim);	copyto!(M,sparsity)
loc2glob = zeros(Int,10);
vertices = zeros(2,3);
J = zeros(2,2);

	
Max_p_tot = zeros(2,2*10^8);
totit=1;

Max_p = zeros(2,1*10^6);
Max_p_new = zeros(2,1*10^6);

	for s = 1:size(Mesh.t_h,2);
	
tid = time();
	
		loc2glob .= Mesh.loc2glob[:,s];
	
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J .= vertices[:,2:end].-vertices[:,1];
	        DetJ = det(J);
	
		Vxᵢ = [V(vertices[:,i]) for i = 1:3];
		diff =  abs(2*Vxᵢ[1]-Vxᵢ[2]-Vxᵢ[3])


		score = 0;

		for i = 1:3 ; for j = 1:2;  if(  abs( mod(vertices[j,i],3) )< 10^-8 ) score+=1; end ;end; end 

		
		if( (diff<eps(10.0)) | (score>=2) )  
			
		add_to_sparse(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)
		
			
		else
	tid_ref = time();	
			Max_p[:,1:3] =vertices;

			no = 1
			
			for i = 1:Nh
					
				#refine
				count = 1;
				new_no = 0;
				for t = 1:no #loop over old, add new
					vrtcs = Max_p[:,(t-1)*3+1:t*3];	
					Vxᵢ = [V(vrtcs[:,i]) for i = 1:3];
					diff =  abs(2*Vxᵢ[1]-Vxᵢ[2]-Vxᵢ[3])
					if(abs(diff)>eps(10.0))
						Max_p_new[:,count:count+12-1] = refine(vrtcs); count+=12; #added 4 simplices 
						new_no +=4;
				
					
					else
						Max_p_new[:,count:count+3-1] = vrtcs; count+=3; 
						new_no+=1;
					end
					
				end
				no = new_no;
				Max_p[:,1:3*no+1] .= Max_p_new[:,1:3*no+1]; 
			end	
			it = 1;
			
			for t in 1:no;
			
				vrtcs = Max_p[:,3*(t-1)+1:3*t];
				Jloc = vrtcs[:,2:3].-vrtcs[:,1];
				dJ = abs(det(Jloc))
			
				New_gp = Jloc*vcat(X_gp',Y_gp');
				New_w = w_gp*dJ;
	
				
				Xr_gp[it:it+N_gp-1] = New_gp[1,:].+vrtcs[1,1]; Yr_gp[it:it+N_gp-1] = New_gp[2,:].+vrtcs[2,1];
				wr_gp[it:it+N_gp-1] = New_w; it+=N_gp;


			end
				add_to_sparse_REF(Xr_gp[1:it-1],Yr_gp[1:it-1],wr_gp[1:it-1],it-1,loc2glob,DetJ,vertices,J)
			
		Max_p_tot[:,totit:totit+3*no-1] = Max_p[:,1:3*no]; totit+=3*no;
			
		end
		
		
		
		
	end
	
	

	
	M+=M'; for i = 1:size(M,1); M[i,i]/=2.; end
	
	return M
end






function Vvᵢvⱼ_piecewise_constant(V,Mesh,Quad,Nh,sparsity)

function add_to_sparse(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)

	M_loc = zeros(10,10);

	for gp = 1:N_gp
		X = X_gp[gp]
		Y = Y_gp[gp]
		w = w_gp[gp]
		x = vertices[:,1]+J*[X;Y];
		
			for i=1:10
				loc_i=vₕ[i](X,Y)
				for j = i:10
					M_loc[j,i]+= loc_i*vₕ[j](X,Y)*V(x)*w*DetJ;
				end
			end 
	end
	
	M_loc+=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
	for i = 1:10
		i_glob = loc2glob[i];
		idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		row = M.rowval[idx];
	
		for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
		
	
end	




function add_to_sparse_REF(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)

	Jinv = 1/det(J)*[ J[2,2] -J[1,2]; -J[2,1] J[1,1]]
	M_loc = zeros(10,10);

	for gp = 1:N_gp
		x = [X_gp[gp]; Y_gp[gp]];
		w = w_gp[gp]
		
		#x = vertices[:,1]+J*[X;Y];
		X = Jinv*(x-vertices[:,1])
			for i=1:10
				loc_i=vₕ[i](X[1],X[2])
				for j = i:10
					M_loc[j,i]+=loc_i*vₕ[j](X[1],X[2])*V(x)*w; 
				end
			end 
	end
	
	M_loc +=M_loc'; for d = 1:10; M_loc[d,d]/=2; end
	for i = 1:10
		i_glob = loc2glob[i];
		idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
		row = M.rowval[idx];
	
		for j = 1:10;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
		m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
		
	
end	



function refine(vertices);
points = [vertices (vertices[:,1]+vertices[:,2])/2 (vertices[:,2]+vertices[:,3])/2 (vertices[:,3]+vertices[:,1])/2];
t = [1,4,6,4,2,5,6,5,3,4,5,6];
return points[:,t];

end



X_gp = Quad["X"];
Y_gp = Quad["Y"];
w_gp = Quad["ω"];

N_gp = length(X_gp);


Xr_gp = zeros(10^6);
Yr_gp = zeros(10^6);
wr_gp = zeros(10^6);

Dim = size(Mesh.p,2);

M = spzeros(Dim,Dim);	copyto!(M,sparsity)
loc2glob = zeros(Int,10);
vertices = zeros(2,3);
J = zeros(2,2);

	
Max_p_tot = zeros(2,2*10^8);
totit=1;

Max_p = zeros(2,1*10^6);
Max_p_new = zeros(2,1*10^6);

	for s = 1:size(Mesh.t_h,2);
	
tid = time();
	
		loc2glob .= Mesh.loc2glob[:,s];
	
		vertices .= Mesh.p[:,Mesh.t_h[:,s]]
	        J .= vertices[:,2:end].-vertices[:,1];
	        DetJ = det(J);
	
		Vxᵢ = [V(vertices[:,i]) for i = 1:3];
		diff =  abs(2*Vxᵢ[1]-Vxᵢ[2]-Vxᵢ[3])



		
		if( (diff<eps(10.0))  )  
			
		add_to_sparse(X_gp,Y_gp,w_gp,N_gp,loc2glob,DetJ,vertices,J)
		
			
		else
	tid_ref = time();	
			Max_p[:,1:3] =vertices;

			no = 1
			
			for i = 1:Nh
					
				#refine
				count = 1;
				new_no = 0;
				for t = 1:no #loop over old, add new
					vrtcs = Max_p[:,(t-1)*3+1:t*3];	
					Vxᵢ = [V(vrtcs[:,i]) for i = 1:3];
					diff =  abs(2*Vxᵢ[1]-Vxᵢ[2]-Vxᵢ[3])
					if(abs(diff)>eps(10.0))
						Max_p_new[:,count:count+12-1] = refine(vrtcs); count+=12; #added 4 simplices 
						new_no +=4;
				
					
					else
						Max_p_new[:,count:count+3-1] = vrtcs; count+=3; 
						new_no+=1;
					end
					
				end
				#println(time()-tid_ref);
				no = new_no;
				#Max_p .= Max_p_new; #ONLY OVERWRITE NEW!
				#overwrite = findfirst(isequal(0),Max_p_new)[2]-1;
				Max_p[:,1:3*no+1] .= Max_p_new[:,1:3*no+1]; 
				#Max_p_new[:,1:overwrite].=0;
			end	
			it = 1;
			
			for t in 1:no;
			
				vrtcs = Max_p[:,3*(t-1)+1:3*t];
				Jloc = vrtcs[:,2:3].-vrtcs[:,1];
				dJ = abs(det(Jloc))
			
				New_gp = Jloc*vcat(X_gp',Y_gp');
				New_w = w_gp*dJ;
	
				
				Xr_gp[it:it+N_gp-1] = New_gp[1,:].+vrtcs[1,1]; Yr_gp[it:it+N_gp-1] = New_gp[2,:].+vrtcs[2,1];
				wr_gp[it:it+N_gp-1] = New_w; it+=N_gp;


			end
				add_to_sparse_REF(Xr_gp[1:it-1],Yr_gp[1:it-1],wr_gp[1:it-1],it-1,loc2glob,DetJ,vertices,J)
			
		#ILLUSTRATE[:,ILLit:ILLit+(it-2)] = [Xr_gp[1:it-1]';Yr_gp[1:it-1]'];ILLit+=it-1;
			#ill = findfirst(Max_p[1,:].==0)-1   ; println(ill," ",totit, " " ,no);
			Max_p_tot[:,totit:totit+3*no-1] = Max_p[:,1:3*no]; totit+=3*no;
			
		end
		
		
		
		
	end
	
	
	#file = matopen("Max_p.mat","w");
	#write(file,"Max_p",Max_p_tot);
	#close(file)
	
	M+=M'; for i = 1:size(M,1); M[i,i]/=2.; end
	
	return M
end









































