module Assemble

export Compute_NL_Energy

using LinearAlgebra,SparseArrays, IterativeSolvers

include("Assemble_Matrices.jl");
include("Assemble_Nonlinear.jl");
include("Assemble_Coarse2Fine.jl");
include("Assemble_SLOD.jl");
include("Basis.jl")
include("Compute_NL_Energy.jl");
include("../../Mesh/3d/3dmesh.jl")
include("../Assemble_vectors_using_tensor.jl");
include("../Assemble_nonlinear_matrix_using_tensor.jl");


function All(mesh,Quad_4,Quad_6,Quad_9,Vs)

	sparsity = MatrixSparsity(mesh)
	
	tid = time()
	A= Assemble.∇vᵢ∇vⱼ(mesh,Quad_4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	Mⱼ = Assemble.vᵢvⱼ(mesh,Quad_6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	tid = time();
	Mv = Assemble.vᵢvⱼ(Vs,mesh,Quad_9,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.P(mesh)
	println("Computed P in ", time()-tid);

#	tid = time()
#	MLz = Assemble.vᵢLzvⱼ(mesh,Quad_6,sparsity)
#	println("Computed Ω in ", time()-tid);
	
	return A ,M,Mv,P,sparsity#,MLz
end



function All(mesh,Quad_4,Quad_6,Quad_9,Vs,Vd);

	sparsity = MatrixSparsity(mesh)

	tid = time();
	∇vᵢ∇vⱼ = Assemble.∇vᵢ∇vⱼ(mesh,Quad_4,sparsity);
	println("Computed stiffness matrix in ", time()-tid);

	tid=time();
	v(x) = 1.0
	vᵢvⱼ = Assemble.vᵢvⱼ(mesh,v,Quad_6,sparsity)
	println("Computed mass matrix in ", time()-tid);

	tid=time();
	#vᵢvⱼVd = Assemble.vᵢvⱼ(mesh,Vd,Quad_9,sparsity);
	vᵢvⱼVd = Assemble.vᵢvⱼ(mesh,Vd,Quad_6,sparsity);
	println("Computed rough matrix in ", time()-tid);

	tid = time();
	vᵢvⱼVs = Assemble.vᵢvⱼ(mesh,Vs,Quad_9,sparsity);
	println("Computed smooth part in ", time()-tid);

	tid = time();
	P = Assemble.Coarse2Fine(mesh)
	println("Computed P in ", time()-tid);

#	tid = time()
#	M_Ω = Assemble.vᵢLzvⱼ(mesh,Quad_6,sparsity)
#	println("Computed Ω in ", time()-tid);

	return ∇vᵢ∇vⱼ ,vᵢvⱼ,vᵢvⱼVs,vᵢvⱼVd,P,sparsity

end



function initial_guess(Ω,mesh,φᵢφⱼ_lu,vᵢvⱼ,dofs_f,V)


	rotation = false;
	if(abs(Ω)>eps(1.0)); rotation = true; end

	if(rotation);
		U = [ (1-Ω)*exp(-dot(Node,Node))+Ω*exp(-dot(Node,Node))*(Node[1]+1im*Node[2]) for Node = eachcol(mesh.p) ]; 
	else 
		U = [ exp(-dot(Node,Node)) for Node = eachcol(mesh.p) ];
	end


	U /= dot(U,vᵢvⱼ*U);
	UGS = φᵢφⱼ_lu\(V* ( (vᵢvⱼ*U)[dofs_f]) );

return UGS;
end



function initial_guess_IterativeSolvers(Ω,mesh,φᵢφⱼ,vᵢvⱼ,dofs_f,V)


	rotation = false;
	if(abs(Ω)>eps(1.0)); rotation = true; end

	if(rotation);
		U = [ (1-Ω)*exp(-dot(Node,Node))+Ω*exp(-dot(Node,Node))*(Node[1]+1im*Node[2]) for Node = eachcol(mesh.p) ]; 
	else 
		U = [ exp(-dot(Node,Node)) for Node = eachcol(mesh.p) ];
	end


	U /= dot(U,vᵢvⱼ*U);
	UGS = cg(φᵢφⱼ,(V* ( (vᵢvⱼ*U)[dofs_f]) ));

return UGS;
end


end
