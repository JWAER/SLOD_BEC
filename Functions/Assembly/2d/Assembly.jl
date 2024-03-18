module Assemble

export Compute_NL_Energy

using LinearAlgebra,SparseArrays

include("Assemble_Matrices.jl");
#include("Assemble_Nonlinear.jl");
include("find_local_bdry.jl")
include("Assemble_SLOD.jl");
include("Basis.jl")
include("Compute_NL_Energy.jl");
include("Assemble_RHS.jl")
include("../../Mesh/2d/GenerateMesh2D_P3.jl")
include("../Assemble_vectors_using_tensor.jl");
include("../Assemble_nonlinear_matrix_using_tensor.jl");


function All(mesh,Quad_4,Quad_6,Quad_9,Vs,Ω)

	sparsity = MatrixSparsity(mesh);
	
	tid = time()
	∇vᵢ∇vⱼ = Assemble.∇vᵢ∇vⱼ(mesh,Quad_4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad_6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	tid = time();
	vᵢvⱼVs = Assemble.Vvᵢvⱼ(Vs,mesh,Quad_9,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.P(mesh)
	println("Computed P in ", time()-tid);

	if(abs(Ω)>eps(1.0));
	tid = time()
	M_Ω = Assemble.M_Ω(mesh,Quad_6,sparsity)
	println("Computed Ω in ", time()-tid);
	else 
	M_Ω = similar(sparsity); fill!(M_Ω,0);
	end
	
	return ∇vᵢ∇vⱼ ,vᵢvⱼ,vᵢvⱼVs,P,M_Ω,sparsity
end






function initial_guess(Ω,mesh,φᵢφⱼ_lu,vᵢvⱼ,dofs_f,V)


	rotation = false;
	if(abs(Ω)>eps(1.0)); rotation = true; end

	if(rotation);
		U = [ ( (1-Ω)*exp(-dot(Node,Node))+Ω*exp(-dot(Node,Node)))*(Node[1]*rand()+1im*Node[2]*rand()) for Node = eachcol(mesh.p) ];
		
		
		# [ (1-Ω)*exp(-dot(Node,Node))+Ω*exp(-dot(Node,Node))*(Node[1]+1im*Node[2]) for Node = eachcol(mesh.p) ]; 
	else 
		U = [ exp(-dot(Node,Node)) for Node = eachcol(mesh.p) ];
	end

#	U /= dot(U,vᵢvⱼ*U);
#	UGS = φᵢφⱼ_lu\(V* ( (vᵢvⱼ*U)[dofs_f] ));


	U = U[dofs_f];
	U /= dot(U,vᵢvⱼ*U);
	UGS = φᵢφⱼ_lu\(V* ( (vᵢvⱼ*U) ));

	
return UGS;
end


end
