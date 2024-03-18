function ϕ(mesh,B,M,P,ℓ)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);

iit = 1;
rowvals = zeros(Int,10^8);
nzval = zeros(10^8);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;


Idx_I = zeros(Int,10^8);
Idx_J = zeros(Int,10^8);
Val = zeros(10^8);




map_C = [findfirst(mesh.CoarseNodes.==dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];


NodalLayer = mesh.CoarseNodalLayer[:,map_C];

	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	#coarse node number map_C[1], is coarse basis 1 in space
	
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];

	
	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
	#setdiff!(Layer,mesh.bdry);

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.Node2Tetrahedra[:,Layer]))),0);
	#fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
	

	#nodes = vec(mesh.loc2glob[:,fine_simplices]);
	#unique!(sort!(nodes));
	
	nodes = vec(mesh.t[:,Coarse_Simplices]);
	unique!(sort!(nodes));
	
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bdry_faces,normals,bdry_tet = find_local_bd(mesh,map_C[C_n]);
	
	bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));

	#CONSTRUCT vᵢ∇ₙvⱼ  = <v_i,dn(v_j)> and all nodes on bdry (not just from coarse scale)

	vᵢ∇ₙvⱼ = Assemble.M_bd(mesh,bdry_faces,normals,bdry_tet,Idx_I,Idx_J,Val);
	
	#remove dofs on bdry to compute "local solutions to P1 input"
	# B_loc = factorize(B[dofs,dofs]);
	B_loc = cholesky(B[dofs,dofs]);
	M_loc = M[dofs,dofs];
	P_loc = Matrix(P[Layer,dofs]);
	N_Coarse = length(Layer);
	φ = [B_loc\(M_loc*P_loc[idx,:]) for idx = 1:N_Coarse];
	
	
	nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMBL DERIVATIVE
	vᵢ∇ₙvⱼ  = vᵢ∇ₙvⱼ[nodes_n,dofs];
  
        vⱼ∇ₙφᵢ=zeros(length(nodes_n),N_Coarse) 
        
	for i = 1:N_Coarse vⱼ∇ₙφᵢ[:,i] = vᵢ∇ₙvⱼ *φ[i]; end;
 
	Least_squares = vⱼ∇ₙφᵢ'*vⱼ∇ₙφᵢ;



	ev, EV = eigen(Least_squares);
	
	idx = findmin(ev)[2];
	if(length(nodes_n)!=length(nodes)); #then contains bdry
		#println(C_n)
		n_idx = findall(Layer.==n)[1];
		#sol  = Phis[n_idx];
		idx = findmax(abs.(EV[n_idx,:]))[2];
		#println(idx, " ", C_n);
	end
	
	g = EV[:,idx]
	sol = sum(φ.*g)


	if(sign(sol[findmax(abs.(sol))[2]]) ==-1) sol*=-1; end
	sol/=maximum(sol);
	
	#ϕ[C_n,dofs] = sol;
	
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol; iit+=len;
	#println(time()-tid)
	end
	
	ϕ = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);
	
	return  sparse(ϕ')
	
	
	
	#return ϕ
end

function find_local_bd(mesh,i)

	nl = setdiff(mesh.CoarseNodalLayer[:,i],0);
	
	tets = setdiff(mesh.Node2Tetrahedra[:,nl],0)
	faces_loc = mesh.t[[6,12,14,15],tets]; #middle is always node no. 6,12,14,15
	faces_loc = vec(faces_loc);
        corresponding_tet = vec([1,1,1,1]*(1:size(tets,1))');
        
        perm = sortperm(faces_loc); faces_loc=faces_loc[perm];
       #bdry = [diff(sort!(faces_loc));1];
        bdry = [diff(faces_loc);1]
        corresponding_tet = corresponding_tet[perm];
        

	it = 1;
	while it<length(bdry)
		
		if(bdry[it]==0); bdry[it+1]=0; it+=2; else it+=1;end
	end

		
	bdry_faces = faces_loc[findall(bdry.>0)];
        bdry_tet = tets[corresponding_tet[findall(bdry.>0)]]
        
        normals = zeros(3,length(bdry_faces));
        
#Compute normals for every face        
for j = 1:length(bdry_faces);
	face = mesh.Faces_dofs[:,bdry_faces[j]];
	n̂ = cross(mesh.p[:,face[4]]-mesh.p[:,face[1]], mesh.p[:,face[10]]-mesh.p[:,face[1]])
	n̂ = n̂/norm(n̂);
	point_on_face = mesh.p[:,face[1]];
	s = 0;
	for i = [1,4,10,20]
	 s+=dot(mesh.p[:,mesh.t[i,bdry_tet[1]]]-point_on_face,n̂);
	end
	if(s>0); n̂ *=-1; end
	normals[:,j] .= n̂ 
end
		
	return bdry_faces, normals,bdry_tet
end
















#------------------------------------------------------ WITH MESH REFINEMENT -----------------------------------------------------------------------





function ϕ_refinement2(mesh,B,M,P,ℓ,sparsity)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);
#ϕ = zeros(length(dofs_C),size(mesh.p,2));
#ϕ = spzeros(length(dofs_C),size(mesh.p,2)); 


#sparse structure of ϕ

iit = 1;
rowvals = zeros(Int,10^8);
nzval = zeros(10^8);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;


#Idx_I = zeros(Int,10^8);
#Idx_J = zeros(Int,10^8);
#Val = zeros(10^8);



map_C = [findfirst(mesh.CoarseNodes.==dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];




NodalLayer = mesh.CoarseNodalLayer[:,map_C];

        for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
        fill!(sparsity.nzval,eps(0.00001));
        
        #println(C_n)
tid = time()
        #coarse node number map_C[1], is coarse basis 1 in space

        n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];


        Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
        #setdiff!(Layer,mesh.bdry);

        Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
        fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end


        nodes = vec(mesh.t[:,fine_simplices]);
        unique!(sort!(nodes));


	setdiff!(Layer,mesh.bdry); #

        bdry_faces,normals,bdry_tet = find_local_bd_refinement(mesh,fine_simplices);

        bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));

        #CONSTRUCT vᵢ∇ₙvⱼ  = <v_i,dn(v_j)> and all nodes on bdry (not just from coarse scale)

	vᵢ∇ₙvⱼ = Assemble.M_bd(mesh,bdry_faces,normals,bdry_tet,sparsity);

        #remove dofs on bdry to compute a-inv. of coarse elements:
        dofs = setdiff(nodes,bdry_nodes_loc)
      #   B_loc = factorize(B[dofs,dofs]);
        B_loc = cholesky(B[dofs,dofs]);
        M_loc = M[dofs,dofs];
        P_loc = Matrix(P[Layer,dofs]);
        N_Coarse = length(Layer);
        φ = [B_loc\(M_loc*P_loc[idx,:]) for idx = 1:N_Coarse];




        nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMBL DERIVATIVE
        vᵢ∇ₙvⱼ  = vᵢ∇ₙvⱼ[nodes_n,dofs];

        vⱼ∇ₙφᵢ=zeros(length(nodes_n),N_Coarse)

        for i = 1:N_Coarse vⱼ∇ₙφᵢ[:,i] = vᵢ∇ₙvⱼ *φ[i]; end;

        Least_squares = vⱼ∇ₙφᵢ'*vⱼ∇ₙφᵢ;



        ev, EV = eigen(Least_squares);

        idx = findmin(ev)[2];
        if(length(nodes_n)!=length(nodes)); #then contains bdry
                #println(C_n)
                n_idx = findall(Layer.==n)[1];
                #sol  = Phis[n_idx];
                idx = findmax(abs.(EV[n_idx,:]))[2];
                #println(idx, " ", C_n);
        end

        g = EV[:,idx]
        sol = sum(φ.*g)


        if(sign(sol[findmax(abs.(sol))[2]]) ==-1) sol*=-1; end
        
        sol /= maximum(sol);


        #ϕ[C_n,dofs] = sol;

        len = length(dofs);
        rowvals[iit:iit+len-1] = dofs#col_s.nzind;
        colptr[C_n+1] = colptr[C_n]+len;
        nzval[iit:iit+len-1] = sol; iit+=len;
   #     println(C_n, " : ",time()-tid)
        end

        ϕ = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);

        return  sparse(ϕ')



        #return ϕ
end











#------------------------------------------------------ Adapted to sum(x.>0) -----------------------------------------------------------------------





function ϕ_refinement(mesh,B,M,P,ℓ,A)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);
#ϕ = zeros(length(dofs_C),size(mesh.p,2));
#ϕ = spzeros(length(dofs_C),size(mesh.p,2)); 


#sparse structure of ϕ

iit = 1;
rowvals = zeros(Int,10^8);
nzval = zeros(10^8);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;


Idx_I = zeros(Int,10^8);
Idx_J = zeros(Int,10^8);
Val = zeros(10^8);




map_C = [findfirst(mesh.CoarseNodes.==dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];


NodalLayer = mesh.CoarseNodalLayer[:,map_C];



#----------------------------Compute a reference canonical SLOD-basis --------------------------------


N =Int( round(length(mesh.CoarseNodes)^(1/3))-2)
if(mod(N,2)==0); origin_C = Int(N^2*N/2-N*N/2+N/2);
else origin_C = Int(ceil(N^3/2));end

n = dofs_C[origin_C] #center node coordinates are mesh.p[:,n];
Layer = setdiff(NodalLayer[:,origin_C],0) #lets do first
Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
nodes = vec(mesh.t[:,fine_simplices]);
unique!(sort!(nodes));
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
bdry_faces,normals,bdry_tet = find_local_bd_refinement(mesh,fine_simplices);
bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));
sol_ref = compute_slod_local(mesh,A,M,P,Layer,n,nodes,bdry_faces,normals,bdry_tet,bdry_nodes_loc,Idx_I,Idx_J,Val);	
println(length(sol_ref));

#--------------------------------------------------------------------------------------------



	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	#println(C_n)
tid = time()	
	#coarse node number map_C[1], is coarse basis 1 in space
	
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];

	
	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
	#setdiff!(Layer,mesh.bdry);

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
	fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
	

	nodes = vec(mesh.t[:,fine_simplices]);
	unique!(sort!(nodes));
	
	#nodes = vec(mesh.t[:,Coarse_Simplices]);
	#unique!(sort!(nodes));
	
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bdry_faces,normals,bdry_tet = find_local_bd_refinement(mesh,fine_simplices);
	
	bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));
        dofs = setdiff(nodes,bdry_nodes_loc)

	println(length(dofs));
	#CONSTRUCT vᵢ∇ₙvⱼ  = <v_i,dn(v_j)> and all nodes on bdry (not just from coarse scale)
	if(   (sum(abs.(mesh.p[:,n] ).< (mesh.H* ℓ))>=1) | (length(dofs)!=length(sol_ref)) );	
	sol = compute_slod_local(mesh,B,M,P,Layer,n,nodes,bdry_faces,normals,bdry_tet,bdry_nodes_loc,Idx_I,Idx_J,Val);	
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol; iit+=len;
	else		
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol_ref; iit+=len;
	end
	
	
	
	end
	
	ϕ = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);

		
	return  sparse(ϕ')
	
	
end

function find_local_bd_refinement(mesh,tets)

	faces_loc = mesh.t[[6,12,14,15],tets]; #middle is always node no. 6,12,14,15
	faces_loc = vec(faces_loc);
        corresponding_tet = vec([1,1,1,1]*(1:size(tets,1))');
        
        perm = sortperm(faces_loc); faces_loc=faces_loc[perm];
       #bdry = [diff(sort!(faces_loc));1];
        bdry = [diff(faces_loc);1]
        corresponding_tet = corresponding_tet[perm];
        

	it = 1;
	while it<length(bdry)
		
		if(bdry[it]==0); bdry[it+1]=0; it+=2; else it+=1;end
	end

		
	bdry_faces = faces_loc[findall(bdry.>0)];
        bdry_tet = tets[corresponding_tet[findall(bdry.>0)]]
        
        normals = zeros(3,length(bdry_faces));
        
#Compute normals for every face        
for j = 1:length(bdry_faces);
	face = mesh.Faces_dofs[:,bdry_faces[j]];
	n̂ = cross(mesh.p[:,face[4]]-mesh.p[:,face[1]], mesh.p[:,face[10]]-mesh.p[:,face[1]])
	n̂ = n̂/norm(n̂);
	point_on_face = mesh.p[:,face[1]];
	s = 0;
	for i = [1,4,10,20]
	 s+=dot(mesh.p[:,mesh.t[i,bdry_tet[1]]]-point_on_face,n̂);
	end
	if(s>0); n̂ *=-1; end
	normals[:,j] .= n̂ 
end
		
	return bdry_faces, normals,bdry_tet
end








#------------------------Compute ϕ CANONICAL-----------------------------------#





function ϕ_canonical(mesh,B,M,P)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);

iit = 1;
rowvals = zeros(Int,10^9);
nzval = zeros(10^9);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;


Idx_I = zeros(Int,10^9);
Idx_J = zeros(Int,10^9);
Val = zeros(10^8);


map_C = [findfirst(mesh.CoarseNodes.==dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];


NodalLayer = mesh.CoarseNodalLayer[:,map_C];


#-------------------------Compute reference basis -----------------------------------
N =Int( round(length(mesh.CoarseNodes)^(1/3))-2)
if(mod(N,2)==0); origin_C = Int(N^2*N/2-N*N/2+N/2);
else origin_C = Int(ceil(N^3/2));end

n = dofs_C[origin_C] #center node coordinates are mesh.p[:,n];
Layer = setdiff(NodalLayer[:,origin_C],0) #lets do first
Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.Node2Tetrahedra[:,Layer]))),0);
nodes = vec(mesh.t[:,Coarse_Simplices]);
unique!(sort!(nodes));
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
bdry_faces,normals,bdry_tet = find_local_bd(mesh,map_C[origin_C]);
bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));
sol_ref = compute_slod_local(mesh,B,M,P,Layer,n,nodes,bdry_faces,normals,bdry_tet,bdry_nodes_loc,Idx_I,Idx_J,Val);	



	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	#coarse node number map_C[1], is coarse basis 1 in space
	
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];

	
	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.Node2Tetrahedra[:,Layer]))),0);
	
	nodes = vec(mesh.t[:,Coarse_Simplices]);
	unique!(sort!(nodes));
	
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bdry_faces,normals,bdry_tet = find_local_bd(mesh,map_C[C_n]);
	
	bdry_nodes_loc = unique(vec(mesh.Faces_dofs[:,bdry_faces]));

        dofs = setdiff(nodes,bdry_nodes_loc)

	nodes_n = setdiff(nodes,mesh.bdry)#
	on_bdry = length(nodes_n)!=length(nodes);
	
	if(on_bdry);	
	sol = compute_slod_local(mesh,B,M,P,Layer,n,nodes,bdry_faces,normals,bdry_tet,bdry_nodes_loc,Idx_I,Idx_J,Val);	
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol; iit+=len;
	else		
	
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol_ref; iit+=len;
	end
	
	
	
	end
	
	ϕ = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);

	
	return  sparse(ϕ')
	
	#return ϕ
end




function compute_slod_local(mesh,B,M,P,Layer,n,nodes,bdry_faces,normals,bdry_tet,bdry_nodes_loc,Idx_I,Idx_J,Val);

	vᵢ∇ₙvⱼ = Assemble.M_bd(mesh,bdry_faces,normals,bdry_tet,Idx_I,Idx_J,Val);
	
	#remove dofs on bdry to compute a-inv. of coarse elements:
	dofs = setdiff(nodes,bdry_nodes_loc)
	# B_loc = factorize(B[dofs,dofs]);
	B_loc = cholesky(B[dofs,dofs]);
	M_loc = M[dofs,dofs];
	P_loc = Matrix(P[Layer,dofs]);
	N_Coarse = length(Layer);
	φ = [B_loc\(M_loc*P_loc[idx,:]) for idx = 1:N_Coarse];
	
	
	nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMBL DERIVATIVE
	vᵢ∇ₙvⱼ  = vᵢ∇ₙvⱼ[nodes_n,dofs];
  
        vⱼ∇ₙφᵢ=zeros(length(nodes_n),N_Coarse) 
        
	for i = 1:N_Coarse vⱼ∇ₙφᵢ[:,i] = vᵢ∇ₙvⱼ *φ[i]; end;
 
	Least_squares = vⱼ∇ₙφᵢ'*vⱼ∇ₙφᵢ;



	ev, EV = eigen(Least_squares);
	
	idx = findmin(ev)[2];
	if(length(nodes_n)!=length(nodes)); #then contains bdry
		#println(C_n)
		n_idx = findall(Layer.==n)[1];
		#sol  = Phis[n_idx];
		idx = findmax(abs.(EV[n_idx,:]))[2];
		#println(idx, " ", C_n);
	end
	
	g = EV[:,idx]
	sol = sum(φ.*g)


	if(sign(sol[findmax(abs.(sol))[2]]) ==-1) sol*=-1; end
	
	sol/= maximum(sol);
	
	return sol
end







