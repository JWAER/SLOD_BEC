function get_idx(x,N)
Nx = Int(round(3N*x[1]+1));
Ny = Int(round(3N*x[2]+1));
Nz = Int(round(3N*x[3]+1));
return Nx+(3N+1)*(Ny-1)+(3N+1)^2*(Nz-1);
end


function generate(N,box_size,ℓ);

p,t,CoarseNodes = generate_P3_mesh(N);

#------------------ Compute Layers of order ℓ -------------------#
Layer = zeros(Int,(2ℓ+1)^3,(N+1)^3);
dx = 1/N;
M = get_idx([1,1,1],N);
for n = 1:length(CoarseNodes);
    count = 1;
    xₙ = p[:,CoarseNodes[n]];

    for nz = -ℓ:ℓ
    	z = xₙ[3]+nz*dx
        for ny = -ℓ:ℓ
            y = xₙ[2]+ny*dx
            for nx = -ℓ:ℓ
              	x = xₙ[1]+nx*dx
		if((x>-eps(0.1))&(abs(x)<1)&(y>-eps(0.1))&(abs(y)<1)&(z>-eps(0.1))&(abs(z)<1))
                idx=  get_idx(xₙ+[nx*dx; ny*dx; nz*dx],N);Layer[count,n]=idx;count+=1;
                end
                #if(Layer[count-1,n]<1 | Layer[count-1,n]>M); count-=1; Layer[count,n] = 0; end 
            end
        end
    end


end
node2tetrahedra = Node2Tetrahedra(t,3N+1);
#---------------Compute list of faces and face to tetrahedra mapping-------
NTet = size(t,2);
Nfaces = 2NTet+N^2*6;
#faces2tet = zeros(Int,2,Nfaces);
#tet2faces = zeros(Int,4,Ntet);
faces_dofs = zeros(Int,10,size(p,2));
faces = zeros(Int,Nfaces*10); #Middle degree of freedom uniquely defines face
it = 1;
for tn = 1:size(t,2);

   faces[it:it+3] = t[[6,12,14,15],tn];it+=4; 
        
end

unique!(sort!(faces))
faces = faces[2:end]; #
face_idx = [6,12,14,15];
local_faces =hcat(1:10, [1,2,3,4,11,12,13,17,18,20],[1,5,8,10,11,14,16,17,19,20], [4,7,9,10,13,15,16,18,19,20])';
for tn = 1:size(t,2);

    for i = 1:4
        face = t[face_idx[i],tn];
        idx = searchsorted(faces,face)[1]
        faces_dofs[:,face] .= t[local_faces[i,:],tn];
    end

end
#------------------------------------------------------------

bdry = zeros(Int,(3N+1)^2*6-(3N+1)*12+8);
it = 1;
for i = 1:size(p,2)
    if( (abs(p[1,i])<eps(1.0)) | (abs(p[1,i]-1)<eps(1.0)) ) 
        bdry[it] = i; it+=1;
    elseif( (abs(p[2,i])<eps(1.0)) | (abs(p[2,i]-1)<eps(1.0)))
        bdry[it] = i ; it+=1;
    elseif( (abs(p[3,i])<eps(1.0)) | (abs(p[3,i]-1)<eps(1.0)))
        bdry[it] = i; it+=1;
    end
end

dofs = setdiff(1:size(p,2),bdry);

H = 2*box_size/N;
p = 2*box_size*(p.-[0.5;0.5;0.5])


return Mesh(H,p,t,box_size,CoarseNodes,Layer,node2tetrahedra,faces,faces_dofs,bdry,dofs)

end












#-------------------------------------With mesh refinement----------------------------------------





function generate(N,Nh,box_size,ℓ);

p,t,not_needed = generate_P3_mesh(N*Nh); #

t_H,CoarseNodes = generate_P1_mesh(N,Nh); #


#CoarseNodes = copy(vec(t_H));
#unique!(sort!(CoarseNodes)); #

p_ref,t_ref = generate_refinement(Nh); #refinement of reference simplex

#------------------ Compute Layers of order ℓ -------------------#
Layer = zeros(Int,(2ℓ+1)^3,(N+1)^3);
dx = 1/N;
M = get_idx([1,1,1],N);
for n = 1:length(CoarseNodes);
    count = 1;
    xₙ = p[:,CoarseNodes[n]];

    for nz = -ℓ:ℓ
    	z = xₙ[3]+nz*dx
        for ny = -ℓ:ℓ
            y = xₙ[2]+ny*dx
            for nx = -ℓ:ℓ
              	x = xₙ[1]+nx*dx
		if((x>-eps(0.1))&(abs(x)<1)&(y>-eps(0.1))&(abs(y)<1)&(z>-eps(0.1))&(abs(z)<1))
                idx=  get_idx(xₙ+[nx*dx; ny*dx; nz*dx],N*Nh);Layer[count,n]=idx;count+=1;
                end
                #if(Layer[count-1,n]<1 | Layer[count-1,n]>M); count-=1; Layer[count,n] = 0; end 
            end
        end
    end


end

FineSimplicesInCoarse = [i:i for i = 1:size(t_H,2)]; #preallocate

#---------Compute all fine tetrahedra using the reference refinement------------
fill!(t,0)
count = 1;
for it_tH = 1:size(t_H,2);
	start = count;
	nodes = p[:,t_H[:,it_tH]];
	p₀ = nodes[:,1];
	J = nodes[:,2:4].-p₀;
	
	p̄ = J*p_ref.+p₀
	for tloc = 1:Nh^3;
		ploc = p̄[:,t_ref[:,tloc]]
		for i = 1:20; t[i,count] = get_idx(ploc[:,i],Nh*N); end
		count +=1;
	end
	FineSimplicesInCoarse[it_tH] = start:count-1;
end
#-------------------------------------------------------------------------------




#FineNode2FineTetrahedra = Node2Tetrahedra(t,3*N*Nh+1); #zeros(Int,24,(3N+1)^3);
CoarseNode2CoarseTetrahedra = Node2Tetrahedra(t_H,3*Nh*N+1);

#---------------Compute list of faces and face to tetrahedra mapping-------
NTet = size(t,2);
Nfaces = 2NTet+(N*Nh)^2*6;
#faces2tet = zeros(Int,2,Nfaces);
#tet2faces = zeros(Int,4,Ntet);
faces_dofs = zeros(Int,10,size(p,2));
faces = zeros(Int,Nfaces*10); #Middle degree of freedom uniquely defines face
it = 1;
for tn = 1:size(t,2);

  #  face1 = t[1:10,tn]; #normal [0,0,-1]
  #  face2 = t[[1,2,3,4,11,12,13,17,18,20],tn]; #normal [0,-1,0]
  #  face3 = t[[4,7,9,10,13,15,16,18,19,20],tn]; #[1,1,1]
  #  face4 = t[[1,5,8,10,11,14,16,17,19,20],tn]; #normal [-1,0,0]

   # face1H = t[[1,4,10],tn];
   # face2H = t[[1,4,20],tn];
   # face3H = t[[4,10,20],tn];
   # face4H = t[[1,10,20],tn];
   faces[it:it+3] = t[[6,12,14,15],tn];it+=4; 
        
end

unique!(sort!(faces))
faces = faces[2:end]; #remove 0
face_idx = [6,12,14,15];
local_faces =hcat(1:10, [1,2,3,4,11,12,13,17,18,20],[1,5,8,10,11,14,16,17,19,20], [4,7,9,10,13,15,16,18,19,20])';
for tn = 1:size(t,2);

    for i = 1:4
        face = t[face_idx[i],tn];
        idx = searchsorted(faces,face)[1]
        faces_dofs[:,face] .= t[local_faces[i,:],tn];
    end

end
#------------------------------------------------------------

bdry = zeros(Int,(3N*Nh+1)^2*6-(3N*Nh+1)*12+8);
it = 1;
for i = 1:size(p,2)
    if( (abs(p[1,i])<eps(1.0)) | (abs(p[1,i]-1)<eps(1.0)) ) 
        bdry[it] = i; it+=1;
    elseif( (abs(p[2,i])<eps(1.0)) | (abs(p[2,i]-1)<eps(1.0)))
        bdry[it] = i ; it+=1;
    elseif( (abs(p[3,i])<eps(1.0)) | (abs(p[3,i]-1)<eps(1.0)))
        bdry[it] = i; it+=1;
    end
end

dofs = setdiff(1:size(p,2),bdry);

H = 2*box_size/N;
p = 2*box_size*(p.-[0.5;0.5;0.5])


return Mesh_refined(H,p,t,t_H,box_size,CoarseNodes,Layer,CoarseNode2CoarseTetrahedra,FineSimplicesInCoarse,faces,faces_dofs,bdry,dofs)

end





















function generate_P3_mesh(N);

	p_ref = Array([0 0 0 
	    1 0 0 
	    0 1 0
	    1 1 0
	    0 0 1
	    1 0 1
	    0 1 1
	    1 1 1 
	    ]');
	t_ref = Array([1 5 6 8
	    1 6 2 8
	    1 2 4 8
	    8 5 7 1
	   1 4 3 8
	   1 3 7 8  ]');

	ref_points = zeros(3,20);
	it =1;
	for zi = 0:3
	    for yi = 0:3
		for xi = 0:3
		    if(xi+yi+zi<=3)
		        ref_points[:,it] = [xi/3;yi/3;zi/3]; 
		        it = it+1;
		    
		    end
		    
		end
	    end
	end

	#N = 2; #will give N^3 boxes
	count = 1;

	t = zeros(Int,20,6*N^3);
	p = zeros(3,(3N+1)^3);

	CoarseNodes = zeros(Int,(N+1)^3);

	dx = 1/N;
	for iz = 1:N+1
	for iy = 1:N+1
	for ix = 1:N+1
	    CoarseNodes[ix+(iy-1)*(N+1)+(iz-1)*(N+1)^2] = get_idx([(ix-1)*dx,(iy-1)*dx,(iz-1)*dx],N);
	end
	end
	end





	dx = 1/(3N)
	for iz = 1:3N+1
	for iy = 1:3N+1
	for ix = 1:3N+1
	p[:,ix+(iy-1)*(3N+1)+(iz-1)*(3N+1)^2] .= [(ix-1)*dx,(iy-1)dx,(iz-1)*dx]; 
	end
	end
	end


	for iz = 1:N
	for iy = 1:N
	for ix = 1:N
		p0 = [1/N*(ix-1);1/N*(iy-1);1/N*(iz-1)]; #translate reference cube here.
		for t_idx = 1:6 # 6 local tetrahedra
			p_loc = p_ref[:,t_ref[:,t_idx]];
			J = hcat(p_loc[:,2]-p_loc[:,1], p_loc[:,3]-p_loc[:,1],p_loc[:,4]-p_loc[:,1]);
			p_t = p0 .+ (p_loc[:,1].+J*ref_points)/N;
			for i = 1:20; t[i,count] = get_idx(p_t[:,i],N);end; count+=1;
		
		end
	end
	end
	end

	return p,t,CoarseNodes


end









function generate_P1_mesh(N,Nh);

	p_ref = Array([0 0 0 
	    1 0 0 
	    0 1 0
	    1 1 0
	    0 0 1
	    1 0 1
	    0 1 1
	    1 1 1 
	    ]');
	t_ref = Array([1 5 6 8
	    1 6 2 8
	    1 2 4 8
	    8 5 7 1
	   1 4 3 8
	   1 3 7 8  ]');


	count = 1;

	t = zeros(Int,4,6*N^3);
	#p = zeros(3,(3N+1)^3);

	CoarseNodes = zeros(Int,(N+1)^3);

	dx = 1/N;
	for iz = 1:N+1
	for iy = 1:N+1
	for ix = 1:N+1
	    CoarseNodes[ix+(iy-1)*(N+1)+(iz-1)*(N+1)^2] = get_idx([(ix-1)*dx,(iy-1)*dx,(iz-1)*dx],N*Nh);
	end
	end
	end




	for iz = 1:N
	for iy = 1:N
	for ix = 1:N
		p0 = [1/N*(ix-1);1/N*(iy-1);1/N*(iz-1)]; #translate reference cube here.
		for t_idx = 1:6 # 6 local tetrahedra
			p_loc = p_ref[:,t_ref[:,t_idx]];
			#J = p_loc[:,2:4].-p_loc[:,1];
			p_t = p0 .+  p_loc/N;
			for i = 1:4; t[i,count] = get_idx(p_t[:,i],N*Nh);end; count+=1;
		
		end
	end
	end
	end

	return t,CoarseNodes


end














function generate_refinement(Nh)

fine_t = zeros(Int,20,Nh^3);

#Generate the refinment of a standard coarse simplex, the rest are just affine mappings of this
	p_ref = Array([0 0 0 
	    1 0 0 
	    0 1 0
	    1 1 0
	    0 0 1
	    1 0 1
	    0 1 1
	    1 1 1 
	    ]');
	t_ref = Array([1 5 6 8
	    1 6 2 8
	    1 2 4 8
	    8 5 7 1
	   1 4 3 8
	   1 3 7 8  ]');


#start with reference discretization

p,t = generate_P3_mesh(Nh);

#J maps reference simplex to simplex number 1
J = hcat(p_ref[:,2],p_ref[:,4],p_ref[:,8]); Jinv = inv(J);

Jinv = [1. -1. 0. ; 0. 1.0 -1. ; 0. 0. 1.0];
count = 1;
for i = 1:size(t,2);
	vertices = p[:,t[ [1,4,10,20],i]];
	
	center = sum(vertices,dims=2)/4;
	X = Jinv*center;
	inside = true;
	
	for i = 1:3
	if( ( X[i]<=-eps(2.) ) | (X[i]>= (1+eps(2.0)) )); inside = false; end
	end
	
	if( sum(X)>=(1.0+eps(2.)) ); inside = false; end
	if(inside); fine_t[:,count].=t[:,i]; count+=1; end

end


#println(size(unique!(vec(fine_t))));
active_nodes = copy(vec(fine_t));
unique!(sort!(active_nodes));

fine_t_ = similar(fine_t); 
for it = 1:length(fine_t); new_idx = searchsorted(active_nodes,fine_t[it])[1]; fine_t_[it] = new_idx; end

p_fine = Jinv*p[:,active_nodes];
return p_fine,fine_t_
end



function Node2Tetrahedra(t,N)
	N2T = zeros(Int,24,N^3);
	#-----------------Compute node2tetrahedra mapping----------------------
	for tn = 1:size(t,2);

	    for i = 1:size(t,1)#20;

		node = t[i,tn];
		current_list = N2T[:,node];
		
		new_list = union(tn,setdiff(current_list,0));
		N2T[1:length(new_list),node] = new_list; 

		end

	end
return N2T
end
