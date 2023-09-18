#Bit of a patchwork. Basically finds the boundary of a square region ... (with cut off corners). 

function find_local_bd(nodes,Mesh,ell);

	#bdry = maximum(Mesh.p[1,:]);
	bdry = Mesh.box;
	x_max = maximum(Mesh.p[1,nodes]);
	x_min = minimum(Mesh.p[1,nodes]);
	
	y_max = maximum(Mesh.p[2,nodes]);
	y_min = minimum(Mesh.p[2,nodes]);
	
	change_row = (3*Mesh.N+1);
	
	l = ell+1; #need one more layer to make it work

	start = nodes[1];#minimum(nodes); 
	count =1;
	
	bd_s = zeros(Int,4*(2*l)+1);

	bd_s[1] = start;
	
	
	for i = 1:2l #ALONG X to the RIGHT
		next_node = bd_s[count]+3*Mesh.Nh; bd_s[count+1] = next_node;count+=1;
		if( abs(Mesh.p[1,next_node]-x_max)<eps(5.) ) break; end
	end
	
	for i = 1:2l #Along y upwards
		next_node = bd_s[count]+3*Mesh.Nh*change_row;bd_s[count+1] = next_node;count+=1;
		if( abs(Mesh.p[2,next_node]-y_max)<eps(5.) )  ;break; end
	end
				
	for i = 1:2l #Along x to the left!
		next_node = bd_s[count]-3*Mesh.Nh;bd_s[count+1] = next_node;count+=1;
		if( abs(Mesh.p[1,next_node]-x_min)<eps(5.) ) break; end
	end
	
	for i = 1:2l #Along y downwards!
		next_node = bd_s[count]-3*Mesh.Nh*change_row;bd_s[count+1] = next_node;count+=1;
		if( abs(Mesh.p[2,next_node]-y_min)<eps(5.) ) break; end
	end
	
	bd_s = intersect(bd_s[1:count-1],nodes);
	bd_fine = zeros(Int,length(bd_s)*Mesh.Nh*3+1);
	idx = 1;
	bd_fine[idx] = bd_s[1];
	for i = 1:length(bd_s); #All degrees of freedom on local bdry
		N_subdiv = Mesh.Nh*3;
		x0 = Mesh.p[:,bd_s[i]]; x1 = Mesh.p[:,bd_s[mod(i,length(bd_s))+1]];
		for k = 1:N_subdiv; coord= (Mesh.box.+(x0+(x1-x0)/N_subdiv*k))/(2*Mesh.box); bd_fine[idx+k] = Coord2Number(coord[1],coord[2],Mesh.N); end; idx += N_subdiv;	
	end

	return bd_fine[1:end-1];#bd_s#Should preserve order

end
