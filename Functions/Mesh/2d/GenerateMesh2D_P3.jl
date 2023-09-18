#Generates a P3-mesh such that dofs correspond to a cartesian grid

#=
*
10
*----*
8    9
*----*----*
5    6    7
*----*----*----*
1    2    3     4
=#
#------------------------------------------------------------------------------
function Coord2Number(x,y,N)
	Nx = round(Int,x*3N)+1;
	Ny = round(Int,y*3N)+1;
	
	return (Ny-1)*(3N+1)+Nx

end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function addN2S(n,s,node2simplex,n_s_per_node)
	for node in n
		if(!(s in  node2simplex[:,node])); node2simplex[n_s_per_node[node],node]=s; n_s_per_node[node]+=1; end
	end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function addloc2glob(x0,x1,x2,n,loc2glob,H,N)
	for i = 1:4; x = x0+(i-1)/3*(x1-x0); loc2glob[i,n] = Coord2Number(x[1],x[2],N); end
  	for i = 5:7; x=x0+(x2-x0)/3+(i-5)/3*(x1-x0); loc2glob[i,n] = Coord2Number(x[1],x[2],N); end
  	for i = 8:9; x=x0+2*(x2-x0)/3+(i-8)/3*(x1-x0); loc2glob[i,n] = Coord2Number(x[1],x[2],N); end
  	loc2glob[10,n] = Coord2Number(x2[1],x2[2],N);
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function Compute_Layers(NH,N,ell);
H = 1/NH;
Layer = zeros((2*ell+1)^2,(NH+1)^2);
CoarseNodes = zeros(Int,(NH+1)^2);
idx_C = 1
	for i_y = 1:NH+1
	for i_x = 1:NH+1
		idx_L = 1;
		x_C = (i_x-1)*H; y_C = (i_y-1)*H; 
		CoarseNodes[idx_C] = Coord2Number(x_C,y_C,N)
		for l_y = -ell:ell#1:(2*ell+1)^2
			y = y_C+l_y*H;
			for l_x = -ell:ell
				x = x_C+l_x*H;
				if(x>-eps(5.) && x<1+eps(5.) && y>-eps(5.) && y<1+eps(5.)) 
					Layer[idx_L,idx_C] = Coord2Number(x,y,N);idx_L+=1;
				end
			end
		end
		idx_C+=1;
	end
	end
return Layer,CoarseNodes
end

#---------------------------------Refine Coarse Simplex --------------------------------------
function Refine_Coarse(Nh,t_h,x0,x1,x2,h_ref,N,count_simplex_fine,loc2glob,node2simplex_fine,n_s_per_node_fine)

#Refine reference simplex then map it to correct place in space using x(X) = x0+(x1-x0)*X+(x2-x0)*Y
h_loc = 1/Nh;
F(X,Y) = x0+(x1-x0)*X+(x2-x0)*Y;
	for j = 0:Nh-1
		Y = j/Nh
		for i = 0:Nh-1-j
			X = i/Nh
			
			idx = count_simplex_fine[1] #mutable object, if not array... 
			#do lower triangle 
			x0_h = F(X,Y)
			x1_h = F(X+h_loc,Y);
			x2_h = F(X,Y+h_loc);
			
			t_h[1,idx] = Coord2Number(x0_h[1],x0_h[2],N); 
			t_h[2,idx] = Coord2Number(x1_h[1],x1_h[2],N);
			t_h[3,idx] = Coord2Number(x2_h[1],x2_h[2],N);
			addloc2glob(x0_h,x1_h,x2_h,idx,loc2glob,h_ref,N);
			addN2S(t_h[:,idx],idx,node2simplex_fine,n_s_per_node_fine)
			count_simplex_fine.+=1;
			
			if(i<Nh-1-j) #DO UPPER
			idx = count_simplex_fine[1] #mutable object, if not array... 
			
			x0_h = F(X+h_loc,Y+h_loc)
			x1_h = F(X,Y+h_loc)
			x2_h = F(X+h_loc,Y);
			
			t_h[1,idx] = Coord2Number(x0_h[1],x0_h[2],N);
			t_h[2,idx] = Coord2Number(x1_h[1],x1_h[2],N);
			t_h[3,idx] = Coord2Number(x2_h[1],x2_h[2],N);
			addloc2glob(x0_h,x1_h,x2_h,idx,loc2glob,h_ref,N);
			addN2S(t_h[:,idx],idx,node2simplex_fine,n_s_per_node_fine)
			
			count_simplex_fine.+=1;
			
			end
		end
	end
end




function Generate(NH,Nh,box,ell);
	#First generate fine mesh on (0,1)^2 then scale and translate
	H_ref = 1/NH;
	h_ref = H_ref/Nh;
	N = NH*Nh;

	#--------  Allocate mesh structures -----------
	t_h = zeros(Int,3,2*N^2);
	FineSimplicesInCoarse= [0:1 for i = 1:2*NH^2]; #initialize Vector{UnitRange{Int64}}
	
	bd = zeros(Int,4*3N);
	loc2glob = zeros(Int,10,2N^2); #Number of simplices is 2N^2
	count_simplex = 1;
	node_count = 1;
	
	node2simplex_fine = zeros(Int,6,(3N+1)^2); #Maps coarse node to coarse simplices containing node
	n_s_per_node_fine = ones(Int, (3N+1)^2); #ONLY FOR COARSE
	
	#----------- Allocate Coarse Mesh --------------------
	t_H = zeros(Int,3,2*NH^2);

	node2simplex_coarse = zeros(Int,6,(3N+1)^2); #Maps coarse node to coarse simplices containing node
	n_s_per_node_coarse = ones(Int, (3N+1)^2); #ONLY FOR COARSE
	

#-------------------Generate nodes (cartesian grid) ----------------------
	p_x = zeros(3N+1,3N+1);
	for i = 1:3N+1
		p_x[:,i] = LinRange(0,1,3N+1);
	end
	p = zeros(2,(3N+1)^2);
	p[1,:] = reshape(p_x,(3N+1)^2);
	p[2,:] = reshape(p_x',(3N+1)^2);
#-----------------------------------------------------------------------#


	
#------------------ Generate coarse simplices ------------------
	for j = 1:NH
		y = (j-1)*H_ref;
		for i = 1:NH
			x = (i-1)*H_ref;
			#do lower triangle 
			t_H[1,count_simplex] = Coord2Number(x,y,N);
			t_H[2,count_simplex] = Coord2Number(x+H_ref,y,N);
			t_H[3,count_simplex] = Coord2Number(x,y+H_ref,N);
			addN2S(t_H[:,count_simplex],count_simplex,node2simplex_coarse,n_s_per_node_coarse); #list of simplices containing node
			count_simplex+=1;
		
			#do upper
			t_H[1,count_simplex] = Coord2Number(x+H_ref,y+H_ref,N);
			t_H[2,count_simplex] = Coord2Number(x,y+H_ref,N);
			t_H[3,count_simplex] = Coord2Number(x+H_ref,y,N);
			addN2S(t_H[:,count_simplex],count_simplex,node2simplex_coarse,n_s_per_node_coarse); #list of simplices containing node
			count_simplex+=1;  
		end
	end


#------------------Generate Fine simplices in accordance to Coarse simplices  ------------------------#
count_simplex_fine = [1];
start = 1;
for s in 1:size(t_H,2);
	x0 = p[:,t_H[1,s]];
	x1 = p[:,t_H[2,s]];
	x2 = p[:,t_H[3,s]];
	Refine_Coarse(Nh,t_h,x0,x1,x2,h_ref,N,count_simplex_fine,loc2glob,node2simplex_fine,n_s_per_node_fine); 
	FineSimplicesInCoarse[s] = start:start+Nh^2-1; start +=Nh^2;
end	
#------------------ DO BOUNDARY, COUNTER CLOCKWISE------------------------------
	for i = 1:4
		bd[1:3N] = 1:3N;
		bd[3N+1:6N] = 3N+1 : 3N+1 : (3N+1)^2-1;	
		bd[6N+1:9N] = (3N+1)^2:-1:(3N+1)^2-3N+1;
		bd[9N+1:12N] = reverse(1+(3N+1):(3N+1):(3N+1)^2-3N)
	end
	
#------------------Compute Layers for each coarse node ----------------------
Layers,CoarseNodes = Compute_Layers(NH,N,ell);
#------------------------- Translate and scale --------------------------------
H = round(10^6*H_ref*box)/10^6;

p = (p.-[0.5;0.5])*box	
	
return Mesh(H,t_H,Nh,N,box/2,t_h,p,loc2glob,CoarseNodes,node2simplex_coarse,node2simplex_fine,Layers,FineSimplicesInCoarse,bd);
end
