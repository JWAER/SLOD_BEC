include("../Functions/Dependencies_3d.jl");
using MKLSparse

β = 100;
Ω = 0.
N  = 8;
ℓ = 1;
box_size = 6; #will give a cube domain from -box to box
ϵ = 1/2;
mesh = Mesh3d.generate(N,box_size,ℓ)
sparsity = Assemble.MatrixSparsity(mesh);

tid_all = time();


@time ∇vᵢ∇vⱼ  = Assemble.∇vᵢ∇vⱼ(mesh,Quad4,sparsity)
v(x) = 1.;
@time vᵢvⱼ = Assemble.vᵢvⱼ(mesh,Quad6,sparsity);

V(x) = 1*dot(x,x)+100*dot(sin.(pi*x),sin.(pi*x));

dynamics = true
Vₜ(x) = 2*dot(x,x);



tid = time();
Vvᵢvⱼ = Assemble.vᵢvⱼ(mesh,V,Quad8,sparsity);
println("computed potential ", time()-tid);
P = Assemble.Coarse2Fine(mesh);




tid = time()
ϕ = Assemble.ϕ_canonical(mesh,∇vᵢ∇vⱼ,vᵢvⱼ,P);
println("computed SLOD-basis ", time()-tid);

ϕ = ϕ[:,mesh.dofs];

∇φᵢ∇φⱼ = ϕ*(∇vᵢ∇vⱼ[mesh.dofs,mesh.dofs]*ϕ')
φᵢφⱼ = ϕ*(vᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ')
Vφᵢφⱼ  = ϕ*(Vvᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');
Vₜvᵢvⱼ = Assemble.vᵢvⱼ(mesh,Vₜ,Quad8,sparsity);
Vₜφᵢφⱼ = ϕ*(Vₜvᵢvⱼ[mesh.dofs,mesh.dofs]*ϕ');
	

∇vᵢ∇vⱼ = 0;
Vvᵢvⱼ = 0; 

SpaceDim_C = size(ϕ,1);


save("./Results/Matrices.jld","ϕ",ϕ,"vᵢvⱼ",vᵢvⱼ,"Vφᵢφⱼ",Vφᵢφⱼ,"Vₜφᵢφⱼ",Vₜφᵢφⱼ,"∇φᵢ∇φⱼ",∇φᵢ∇φⱼ,"φᵢφⱼ",φᵢφⱼ)

