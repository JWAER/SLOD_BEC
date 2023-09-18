function MatrixSparsity(mesh)
	N = size(mesh.p,2);
	M = spzeros(N,N);
	
	Idx_I = ones(Int,10^7);
   	Idx_J = ones(Int,10^7);
   	idx_IJ = 1;
   	
   	loc2glob = zeros(Int,20);
   
for t = 1:size(mesh.t,2);


        loc2glob .= mesh.t[:,t];

	    for i = 1:20
	       i_glob = loc2glob[i];
	        for j = 1:20;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
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




















function  vᵢvⱼ(  mesh,weight,Quad,sparsity)
#Calculates <weight[x]* phi_j, phi_i>
N_gp = Int(Quad["N_gp"][1])

N = size(mesh.p,2);
   #Preallocate
   M = spzeros(N,N); copyto!(M,sparsity);
   Vertices = zeros(3,4);
   J = zeros(3,3);
   loc2glob = zeros(Int,20);
   simplex = zeros(Int,4);
   M_loc = zeros(20,20);
   
for t = 1:size(mesh.t,2);

	fill!(M_loc,0);
        simplex .= mesh.t[[1,4,10,20],t];

	Vertices .= mesh.p[:,simplex];
        J .= Vertices[:,2:4].-Vertices[:,1];
        #hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1], Vertices[:,4]-Vertices[:,1]);


        DetJ = det(J);
        #DetJ = J[1,1]*J[2,2]*J[3,3]+J[1,2]*J[2,3]*J[3,1]+J[1,3]*J[2,1]*J[3,2]-J[1,1]*J[2,3]*J[3,2]-J[1,2]*J[2,1]*J[3,3]-J[1,3]*J[2,2]*J[3,1];
      
        #---Calculate local matrix-------
        loc2glob .= mesh.t[:,t];

   	for gp = 1:N_gp
	    X = Quad["X"][gp]; Y = Quad["Y"][gp];Z = Quad["Z"][gp];ω =Quad["ω"][gp];
		
	    w = weight(J*[X;Y;Z]+Vertices[:,1])*ω*DetJ;
		  
	    X̄ = [X,Y,Z];
		   
	    for i = 1:20
	    	for j = i:20;
		M_loc[j,i]+=w*vₕ[i](X̄)*vₕ[j](X̄);#M_loc[i,j];
				
		end
	    end
	
	end
	
	M_loc +=M_loc'; for d = 1:20; M_loc[d,d]/=2; end
	
    for i = 1:20
	       i_glob = loc2glob[i];
	       idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
	       row = M.rowval[idx];
		
	        for j = 1:20;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
	        m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
	
	

  
end
M+=M'; for i = 1:size(M,1); M[i,i]/=2;end
  
    return M;
end





function  vᵢvⱼ(  mesh,Quad,sparsity)
#Calculates <weight[x]* phi_j, phi_i>
N_gp = Int(Quad["N_gp"][1])

N = size(mesh.p,2);
   #Preallocate
   M = spzeros(N,N); copyto!(M,sparsity);
   Vertices = zeros(3,4);
   J = zeros(3,3);
   loc2glob = zeros(Int,20);
   simplex = zeros(Int,4);
   M_loc = zeros(20,20);
   
for modular = 1:6
for t = modular:6:size(mesh.t,2);


 simplex .= mesh.t[[1,4,10,20],t];
 loc2glob .= mesh.t[:,t];


if(t == modular)
	fill!(M_loc,0);
       
	Vertices .= mesh.p[:,simplex];
        J .= Vertices[:,2:4].-Vertices[:,1];
        #hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1], Vertices[:,4]-Vertices[:,1]);


        DetJ = det(J);
        #DetJ = J[1,1]*J[2,2]*J[3,3]+J[1,2]*J[2,3]*J[3,1]+J[1,3]*J[2,1]*J[3,2]-J[1,1]*J[2,3]*J[3,2]-J[1,2]*J[2,1]*J[3,3]-J[1,3]*J[2,2]*J[3,1];
      
        #---Calculate local matrix-------
       
   	for gp = 1:N_gp
	    X = Quad["X"][gp]; Y = Quad["Y"][gp];Z = Quad["Z"][gp];ω =Quad["ω"][gp];
		
	    w = ω*DetJ;
		  
	    X̄ = [X,Y,Z];
		   
	    for i = 1:20
	    	for j = i:20;
		M_loc[j,i]+=w*vₕ[i](X̄)*vₕ[j](X̄);#M_loc[i,j];
				
		end
	    end
	
	end
	
	M_loc +=M_loc'; for d = 1:20; M_loc[d,d]/=2; end
end
    for i = 1:20
	       i_glob = loc2glob[i];
	       idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
	       row = M.rowval[idx];
		
	        for j = 1:20;   j_glob = loc2glob[j]; if(j_glob>=i_glob);
	        m = searchsorted(row,j_glob)[1];
		M.nzval[idx[m]]+=M_loc[j,i];
				
	       end end
	end
	
	

end  
end
M+=M'; for i = 1:size(M,1); M[i,i]/=2;end
  
    return M;
end










function ∇vᵢ∇vⱼ_OLD(mesh,Quad,sparsity)
N = size(mesh.p,2);
N_gp = Int(Quad["N_gp"][1])

A = spzeros(N,N);copyto!(A,sparsity);
A_loc = zeros(20,20);
#Calculates du,dv
   #preallocate
   Vertices = zeros(3,4);
   J = zeros(3,3);
   JTinv = zeros(3,3);
   loc2glob = zeros(Int,20);
   simplex = zeros(Int,4);
     ∇φi  = zeros(3);
     ∇φj = zeros(3); 
    X̄ = zeros(3);
    
  #  row = zeros(Int,200);

for t = 1:size(mesh.t,2);


#tid = time();

fill!(A_loc,0.0)

        simplex .= mesh.t[[1,4,10,20],t];

	Vertices .= mesh.p[:,simplex];
        


      J .= Vertices[:,2:4].-Vertices[:,1];
      DetJ = J[1,1]*J[2,2]*J[3,3]+J[1,2]*J[2,3]*J[3,1]+J[1,3]*J[2,1]*J[3,2]-J[1,1]*J[2,3]*J[3,2]-J[1,2]*J[2,1]*J[3,3]-J[1,3]*J[2,2]*J[3,1];
      
      
      JTinv .= inv(J)';	 #

         #Add to Global Matrix
	
        loc2glob .= mesh.t[:,t];
#tid_comp = time();	
		for gp = 1:N_gp
		    X̄[1] = Quad["X"][gp]; X̄[2] = Quad["Y"][gp];X̄[3] = Quad["Z"][gp];ω =Quad["ω"][gp];

		   
		    for i = 1:20
		    	∇φi .= (JTinv*∇vₕ[i](X̄))*ω*DetJ#*DetJ;
		
		        for j = i:20;
		    	    ∇φj .= JTinv*∇vₕ[j](X̄)
		           A_loc[j,i] +=  dot(∇φi,∇φj); 
		         
		        end 
		    end
	

		end
		
		A_loc+=A_loc'; for i = 1:20; A_loc[i,i]/=2; end
	
#println("tid comp ", time()-tid);
#tid = time();
#tid_add = time();

		 for i = 1:20 ;  i_glob = loc2glob[i];
		     if(i_glob!=0)
	            idx = A.colptr[i_glob]:A.colptr[i_glob+1]-1;
	 	    row = A.rowval[idx]; #memory allocation
	 	    #row[1:(idx.stop-idx.start+1)].=A.rowval[idx]
	 	    for j = 1:20; j_glob = loc2glob[j];if(j_glob>= i_glob)
		        m = searchsorted(row,j_glob)[1];#findfirst(isequal(j_glob),row);
		   	A.nzval[idx[m]]+=A_loc[j,i];
		   	#A[j_glob,i_glob] += A_loc[j,i]
		    end end  
         	end end
#println(time()-tid);
    end

A+=A'; for i = 1:size(A,1); A[i,i]/=2;end

        return A;

end





function ∇vᵢ∇vⱼ(mesh,Quad,sparsity)
N = size(mesh.p,2);
N_gp = Int(Quad["N_gp"][1])

A = spzeros(N,N);copyto!(A,sparsity);
A_loc = zeros(20,20);
#Calculates du,dv
   #preallocate
   Vertices = zeros(3,4);
   J = zeros(3,3);
   JTinv = zeros(3,3);
   loc2glob = zeros(Int,20);
   simplex = zeros(Int,4);
     ∇φi  = zeros(3);
     ∇φj = zeros(3); 
    X̄ = zeros(3);
    
  #  row = zeros(Int,200);

for modular = 1:6

for t = modular:6:size(mesh.t,2);


#tid = time();
        simplex .= mesh.t[[1,4,10,20],t];
        loc2glob .= mesh.t[:,t];



if(t == modular) #Compute reference
	fill!(A_loc,0.0);
	Vertices .= mesh.p[:,simplex];
        


      J .= Vertices[:,2:4].-Vertices[:,1];
      DetJ = J[1,1]*J[2,2]*J[3,3]+J[1,2]*J[2,3]*J[3,1]+J[1,3]*J[2,1]*J[3,2]-J[1,1]*J[2,3]*J[3,2]-J[1,2]*J[2,1]*J[3,3]-J[1,3]*J[2,2]*J[3,1];
      
      
      JTinv .= inv(J)';	 #

         #Add to Global Matrix
	
   
#tid_comp = time();	
		for gp = 1:N_gp
		    X̄[1] = Quad["X"][gp]; X̄[2] = Quad["Y"][gp];X̄[3] = Quad["Z"][gp];ω =Quad["ω"][gp];

		   
		    for i = 1:20
		    	∇φi .= (JTinv*∇vₕ[i](X̄))*ω*DetJ#*DetJ;
		
		        for j = i:20;
		    	    ∇φj .= JTinv*∇vₕ[j](X̄)
		           A_loc[j,i] +=  dot(∇φi,∇φj); 
		         
		        end 
		    end
	

		end
		
		A_loc+=A_loc'; for i = 1:20; A_loc[i,i]/=2; end
end
#println("tid comp ", time()-tid);
#tid = time();
#tid_add = time();

		 for i = 1:20 ;  i_glob = loc2glob[i];
		     if(i_glob!=0)
	            idx = A.colptr[i_glob]:A.colptr[i_glob+1]-1;
	 	    row = A.rowval[idx]; #memory allocation
	 	    #row[1:(idx.stop-idx.start+1)].=A.rowval[idx]
	 	    for j = 1:20; j_glob = loc2glob[j];if(j_glob>= i_glob)
		        m = searchsorted(row,j_glob)[1];#findfirst(isequal(j_glob),row);
		   	A.nzval[idx[m]]+=A_loc[j,i];
		   	#A[j_glob,i_glob] += A_loc[j,i]
		    end end  
         	end end
#println(time()-tid);
    end
    
end    

A+=A'; for i = 1:size(A,1); A[i,i]/=2;end

        return A;

end












function  M_bd( Mesh,bdry_faces,normals,bdry_tet,Idx_I,Idx_J,V)
#Calculates

#----------Quad 2D---------------
Quad = Dict{String,Array{Float64}{1}}(
"X" =>  [ 0.33333333333333337
 0.470142064105115
 0.470142064105115
 0.05971587178977
 0.10128650732346
 0.10128650732346
 0.797426985353087],
"Y" =>  [ 0.33333333333333337
 0.470142064105115
 0.05971587178977
 0.470142064105115
 0.10128650732346
 0.797426985353087
 0.10128650732346],
"ω" =>[ 0.1125
 0.066197076394253
 0.066197076394253
 0.066197076394253
 0.0629695902724135
 0.0629695902724135
 0.0629695902724135])
N_gp = 7#Int(Quad["N_gp"][1])
#--------------------------------
N = size(Mesh.p,2);

M = spzeros(N,N);
  # Idx_I = ones(Int,10^5);
  # Idx_J = ones(Int,10^5);
  # V = zeros(10^5);
   idx_IJ = 1;
         p₀ = zeros(3);


        possible_face_dofs_loc = [1 2 3 4 5 6 7 8 9 10;
                         1 2 3 4 11 12 13 17 18 20;
                         1 5 8 10 11 14 16 17 19  20;
                         4 7  9 10 13 15 16 18 19 20]
       possible_Face_Jacobian = [ [1 0 ; 0 1; 0 0 ] , [-1 -1; 0 1 ; 1 0  ] ,[1 0; 0 0 ; 0 1], [0 0 ; 1 0 ; 0 1]]
       


 for it = 1:length(bdry_faces);
        face = bdry_faces[it];
        n̂ = normals[:,it];
        tet = Mesh.t[:,bdry_tet[it]];


        face_simplex = Mesh.Faces_dofs[:,it];

        Vertices = Mesh.p[:,tet[[1,4,10,20]]];


        J = Vertices[:,2:end].-Vertices[:,1];
        #hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1],Vertices[:,4]-Vertices[:,1]);
        JTinv = inv(J)';

        if(face == tet[6]) face_no = 1; elseif(face==tet[12]) face_no = 2; elseif(face==tet[14]); face_no = 3; elseif(face == tet[15]); face_no = 4; end

	face_dofs_loc = possible_face_dofs_loc[face_no,:];
	Face_Jacobian = possible_Face_Jacobian[face_no]
  
  	fill!(p₀,0.)
        if(face_no == 2); p₀.= [1;0;0]; end

        DetJ = abs(det(J));
      
        #---Calculate local matrix-------
      #  loc2glob = tet;
		for gp = 1:N_gp
		    X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

                    X̄ = p₀+Face_Jacobian*[X;Y];
		
		   # w = weight(J*[X;Y;Z]+Vertices[:,1])
		  
		   
		    for j = 1:20 #face_dofs_loc
                        ∇jn = dot(JTinv*∇vₕ[j](X̄),n̂)*ω*DetJ;

		        for i = face_dofs_loc
                            Idx_I[idx_IJ] = tet[i]; Idx_J[idx_IJ] = tet[j]; V[idx_IJ] = ∇jn*vₕ[i](X̄) ; idx_IJ += 1;
		        end

		    end
	

		end

     #---Local matrix calculated------

     #Added   
	if(idx_IJ>0.8*10^8) M+= sparse(Idx_I[1:idx_IJ-1],Idx_J[1:idx_IJ-1],V[1:idx_IJ-1],N,N); fill!(Idx_I,1);fill!(Idx_J,1);fill!(V,0.0); idx_IJ = 1; end
	

end

M+= sparse(Idx_I[1:idx_IJ-1],Idx_J[1:idx_IJ-1],V[1:idx_IJ-1],N,N); 


    return M;
end







function  M_bd( Mesh,bdry_faces,normals,bdry_tet,M)
#Calculates

#----------Quad 2D---------------
Quad = Dict{String,Array{Float64}{1}}(
"X" =>  [ 0.33333333333333337
 0.470142064105115
 0.470142064105115
 0.05971587178977
 0.10128650732346
 0.10128650732346
 0.797426985353087],
"Y" =>  [ 0.33333333333333337
 0.470142064105115
 0.05971587178977
 0.470142064105115
 0.10128650732346
 0.797426985353087
 0.10128650732346],
"ω" =>[ 0.1125
 0.066197076394253
 0.066197076394253
 0.066197076394253
 0.0629695902724135
 0.0629695902724135
 0.0629695902724135])
N_gp = 7#Int(Quad["N_gp"][1])
#--------------------------------
N = size(Mesh.p,2);


         p₀ = zeros(3);


        possible_face_dofs_loc = [1 2 3 4 5 6 7 8 9 10;
                         1 2 3 4 11 12 13 17 18 20;
                         1 5 8 10 11 14 16 17 19  20;
                         4 7  9 10 13 15 16 18 19 20]
       possible_Face_Jacobian = [ [1 0 ; 0 1; 0 0 ] , [-1 -1; 0 1 ; 1 0  ] ,[1 0; 0 0 ; 0 1], [0 0 ; 1 0 ; 0 1]]
       
       
       local_vals = zeros(20,20);
       X̄ = zeros(3);
       


 for it = 1:length(bdry_faces);
        face = bdry_faces[it];
        n̂ = normals[:,it];
        tet = Mesh.t[:,bdry_tet[it]];


        face_simplex = Mesh.Faces_dofs[:,it];

        Vertices = Mesh.p[:,tet[[1,4,10,20]]];


        J = Vertices[:,2:end].-Vertices[:,1];
        #hcat(Vertices[:,2]-Vertices[:,1], Vertices[:,3]-Vertices[:,1],Vertices[:,4]-Vertices[:,1]);
        JTinv = inv(J)';

        if(face == tet[6]) face_no = 1; elseif(face==tet[12]) face_no = 2; elseif(face==tet[14]); face_no = 3; elseif(face == tet[15]); face_no = 4; end

	face_dofs_loc = possible_face_dofs_loc[face_no,:];
	Face_Jacobian = possible_Face_Jacobian[face_no]
  
  	fill!(p₀,0.)
        if(face_no == 2); p₀.= [1;0;0]; end

        DetJ = abs(det(J));
      
        #---Calculate local matrix-------
      #  loc2glob = tet;
   		for gp = 1:N_gp
		    X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

                    X̄ .= p₀+Face_Jacobian*[X;Y];
		
		   # w = weight(J*[X;Y;Z]+Vertices[:,1])
		  
		   
		    for j = 1:20 #face_dofs_loc
                        ∇jn = dot(JTinv*∇vₕ[j](X̄),n̂)*ω*DetJ;

		        for i = face_dofs_loc
		            local_vals[i,j] = ∇jn*vₕ[i](X̄);
                           # Idx_I[idx_IJ] = tet[i]; Idx_J[idx_IJ] = tet[j]; V[idx_IJ] = ∇jn*vₕ[i](X̄) ; idx_IJ += 1;
		        end

		    end
	

		end
     #---Local matrix calculated------
     #=
 	 for i = face_dofs_loc
	       i_glob = tet[i];
	       idx = M.colptr[i_glob]:M.colptr[i_glob+1]-1;
	       row = M.rowval[idx];
		
	        for j = 1:20;  
	        	m = searchsorted(row,tet[j])[1];
			M.nzval[idx[m]]+=local_vals[i,j];
				
	       end end
	end
	=#
 	 for j = 1:20
	       j_glob = tet[j];
	       idx = M.colptr[j_glob]:M.colptr[j_glob+1]-1;
	       row = M.rowval[idx];
		
	        for i = face_dofs_loc;  
	        	m = searchsorted(row,tet[i])[1];
			M.nzval[idx[m]]+=local_vals[i,j];
				
	       end 
	end
	
	fill!(local_vals,0.0);


	

end



    return M;
end





function vᵢLzvⱼ(Mesh,Quad,sparsity)
	Dim = size(Mesh.p,2);
	N_gp = Int(Quad["N_gp"][1]);

	I_vec = zeros(Int,10^6);
	J_vec = zeros(Int,10^6);
	Val_vec = 0im*zeros(10^6);
	idx = 1;
	
	M_Ω = 0im*spzeros(Dim,Dim);	


	for s = 1:size(Mesh.t,2);
		loc2glob = Mesh.t[:,s];
		vertices = Mesh.p[:,Mesh.t[[1,4,10,20],s]]
	        J = hcat(vertices[:,2]-vertices[:,1], vertices[:,3]-vertices[:,1],vertices[:,4]-vertices[:,1]);
		DetJ = det(J);
		JTinv = inv(J)';	
		
		for gp = 1:N_gp
		X = Quad["X"][gp];
		Y = Quad["Y"][gp];
		Z = Quad["Z"][gp];
		w = Quad["ω"][gp]
		yx =  (vertices[:,1]+J*[X;Y;Z]) ;
		yx[1:2].=reverse(yx[1:2]);
		yx[2]*=-1;
		yx[3] = 0. # will give (y, -x, 0)		
			for i=1:20
				grad_i = JTinv*∇vₕ[i]([X,Y,Z]);
				for j = 1:20
				#grad_j = JTinv*∇vₕ[j](X,Y);
					I_vec[idx] = loc2glob[j]; J_vec[idx] = loc2glob[i]; Val_vec[idx] = -dot(yx,grad_i)*vₕ[j]([X,Y,Z])*w*DetJ; idx+=1; #provided dirichlet bcs
				end
			end 
		end
		
		
		if(idx>0.8*10^6); M_Ω+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],Val_vec[1:idx-1],Dim,Dim); idx = 1; end
		
	end
	
	M_Ω+=sparse(I_vec[1:idx-1],J_vec[1:idx-1],Val_vec[1:idx-1],Dim,Dim);
	
#	M_Ω+=M_Ω'; for i = 1:size(M_Ω,1); M_Ω[i,i]/=2.; end
	return M_Ω
end


