function ϕ(mesh,B,M,P,ℓ)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);


#sparse structure of ϕ

iit = 1;
rowvals = zeros(Int,1*10^9);
nzval = zeros(10^9);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;





map_C = [searchsorted(mesh.CoarseNodes,dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];


NodalLayer = mesh.NodalLayer[:,map_C];

	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];

	#origin = intersect(findall(mesh.p[1,:].==0),findall(mesh.p[2,:].==0))[1];

	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
	#setdiff!(Layer,mesh.bdry);

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
	fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
	

	nodes = vec(mesh.loc2glob[:,fine_simplices]);
	unique!(sort!(nodes));
	
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bd = find_local_bd(nodes,mesh,ℓ);

	#CONSTRUCT vᵢ∇ₙvⱼ  = <vᵢ,∇ₙ(vⱼ)> and all nodes on bdry (not just from coarse scale)
	
	vᵢ∇ₙvⱼ = Assemble.M_bc(mesh,bd,bd[1:3:end-1]);
	#remove dofs on bdry to compute a-inv. of coarse elements:
	dofs = setdiff(nodes,bd)
	B_loc = factorize(B[dofs,dofs]);
	#B_loc = cholesky(B[dofs,dofs]);
	M_loc = M[dofs,dofs];
	P_loc = Matrix(P[Layer,dofs]);
	N_Coarse = length(Layer);
	φ = [B_loc\(M_loc*P_loc[idx,:]) for idx = 1:N_Coarse];
	
	
	nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMAL DERIVATIVE
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
	
	len = length(dofs);
	rowvals[iit:iit+len-1] = dofs#col_s.nzind;
	colptr[C_n+1] = colptr[C_n]+len;
	nzval[iit:iit+len-1] = sol; iit+=len;
	
	end
	
	vod = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);
	
	return  sparse(vod')
	
	
	
	#return VLOD
end







function ϕ_canonical(mesh,B,M,P,ℓ)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);

map_C = [searchsorted(mesh.CoarseNodes,dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];
NodalLayer = mesh.NodalLayer[:,map_C];


iit = 1;
rowvals = zeros(Int,10^8);
nzval = zeros(10^8);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;
	
#---------------------------Compute reference basis ---------------------------------------------#
#origin_C = intersect(findall(mesh.p[1,dofs_C].==0),findall(mesh.p[2,dofs_C].==0))[1];
Nn =Int( round(length(mesh.CoarseNodes)^(1/2))-2)
if(mod(Nn,2)==0); origin_C = Int(Nn*Nn/2-Nn/2);
else origin_C = Int(ceil(Nn^2/2));end


n = dofs_C[origin_C] #center node coordinates are mesh.p[:,n];

Layer = setdiff(NodalLayer[:,origin_C],0) #lets do first
Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
nodes = vec(mesh.loc2glob[:,fine_simplices]);
unique!(sort!(nodes));
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!

bd = find_local_bd(nodes,mesh,ℓ);
dofs = setdiff(nodes,bd)
nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMAL DERIVATIVE

sol_ref = compute_slod_local(mesh,B,M,P,bd,Layer,nodes,nodes_n,n,dofs);
col_s = sparse(sol_ref);
		
	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	tid = time();
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];

	#origin = intersect(findall(mesh.p[1,:].==0),findall(mesh.p[2,:].==0))[1];

	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
	#setdiff!(Layer,mesh.bdry);

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
	fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
	

	nodes = vec(mesh.loc2glob[:,fine_simplices]);
	unique!(sort!(nodes));
	
	setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bd = find_local_bd(nodes,mesh,ℓ);
	dofs = setdiff(nodes,bd)
	
	#CONSTRUCT vᵢ∇ₙvⱼ  = <v_i,dn(v_j)> and all nodes on bdry (not just from coarse scale)
	nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMAL DERIVATIVE
	
	on_bdry = length(nodes_n)!=length(nodes);
	
	if(on_bdry);  #then contains bdry
		sol = compute_slod_local(mesh,B,M,P,bd,Layer,nodes,nodes_n,n,dofs);
		#col_s = sparse(sol);
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
	
	
	vod = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);
	
	return  sparse(vod')
end





function compute_slod_local(mesh,B,M,P,bd,Layer,nodes,nodes_n,n,dofs);

	vᵢ∇ₙvⱼ  = Assemble.M_bc(mesh,bd,bd[1:3:end-1])[nodes_n,dofs];
 
	#remove dofs on bdry to compute a-inv. of coarse elements:
	dofs = setdiff(nodes,bd)
	B_loc = factorize(B[dofs,dofs]);

	M_loc = M[dofs,dofs];
	P_loc = Matrix(P[Layer,dofs]);
	N_Coarse = length(Layer);
	φ = [B_loc\(M_loc*P_loc[idx,:]) for idx = 1:N_Coarse];
	
	vⱼ∇ₙφᵢ = zeros(length(nodes_n),N_Coarse) 
	for i = 1:N_Coarse vⱼ∇ₙφᵢ[:,i] = vᵢ∇ₙvⱼ *φ[i]; end;
	 
	Least_squares = vⱼ∇ₙφᵢ'*vⱼ∇ₙφᵢ;

	ev, EV = eigen(Least_squares);
	
	idx = findmin(abs.(ev))[2];
	
	if(length(nodes_n)!=length(nodes)); #then contains bdry
		n_idx = findall(Layer.==n)[1];
		idx = findmax(abs.(EV[n_idx,:]))[2];
	end
	
	g = EV[:,idx]#EV[:,idx];
	sol = sum(φ.*g)


	
	if(sign(sol[findmax(abs.(sol))[2]]) ==-1) sol*=-1; end
	
	sol /= maximum(sol);
	
	return sol;
	

end




































































function ϕ2(mesh,B,M,P,ℓ,A)

dofs_C = setdiff(mesh.CoarseNodes,mesh.bdry);

map_C = [searchsorted(mesh.CoarseNodes,dofs_C[i])[1] for i = 1:length(dofs_C)]; #Map_C i such that dofs_C = mesh.CoarseNodes[map_C];
NodalLayer = mesh.NodalLayer[:,map_C];


iit = 1;
rowvals = zeros(Int,10^8);
nzval = zeros(10^8);
colptr = zeros(Int,length(dofs_C)+1);
colptr[1]=1;
	
#---------------------------Compute reference basis ---------------------------------------------#
#origin_C = intersect(findall(mesh.p[1,dofs_C].==0),findall(mesh.p[2,dofs_C].==0))[1];
Nn =Int( round(length(mesh.CoarseNodes)^(1/2))-2)
if(mod(Nn,2)==0); origin_C = Int(Nn*Nn/2-Nn/2);
else origin_C = Int(ceil(Nn^2/2));end


n = dofs_C[origin_C] #center node coordinates are mesh.p[:,n];

Layer = setdiff(NodalLayer[:,origin_C],0) #lets do first
Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
nodes = vec(mesh.loc2glob[:,fine_simplices]);
unique!(sort!(nodes));
setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!

bd = find_local_bd(nodes,mesh,ℓ);
dofs = setdiff(nodes,bd)
nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMAL DERIVATIVE

sol_ref = compute_slod_local(mesh,A,M,P,bd,Layer,nodes,nodes_n,n,dofs);
col_s = sparse(sol_ref);
		
	for C_n = 1:length(dofs_C)#length(mesh.CoarseNodes);
	tid = time();
	n = dofs_C[C_n] #center node coordinates are mesh.p[:,n];
	x = mesh.p[:,n];

	#origin = intersect(findall(mesh.p[1,:].==0),findall(mesh.p[2,:].==0))[1];

	Layer = setdiff(NodalLayer[:,C_n],0) #lets do first
	#setdiff!(Layer,mesh.bdry);

	Coarse_Simplices=setdiff(unique!(sort!(vec(mesh.CoarseNode2CoarseSimplices[:,Layer]))),0);
	fine_simplices = []; for C_S in Coarse_Simplices; union!(fine_simplices,mesh.FineSimplicesInCoarse[C_S]); end
	

	nodes = vec(mesh.loc2glob[:,fine_simplices]);
	unique!(sort!(nodes));
	
	setdiff!(Layer,mesh.bdry); #need to be here, otherwise some simplices would be omitted!
	
	bd = find_local_bd(nodes,mesh,ℓ);
	dofs = setdiff(nodes,bd)
	
	#CONSTRUCT vᵢ∇ₙvⱼ  = <v_i,dn(v_j)> and all nodes on bdry (not just from coarse scale)
	nodes_n = setdiff(nodes,mesh.bdry)#WHERE TO ENFORCE NORMAL DERIVATIVE
	
	on_bdry = length(nodes_n)!=length(nodes);
	
	if(on_bdry |  (sum(abs.(x).<((ℓ+4)*mesh.H))>0) );  
		sol = compute_slod_local(mesh,B,M,P,bd,Layer,nodes,nodes_n,n,dofs);
		#col_s = sparse(sol);
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
	
	
	vod = SparseMatrixCSC(size(mesh.p,2),length(dofs_C),colptr,rowvals[1:iit-1],nzval[1:iit-1]);
	
	return  sparse(vod')
end



