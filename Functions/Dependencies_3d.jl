using LinearAlgebra, SparseArrays, IterativeSolvers,JLD, MKLSparse, Preconditioners


include("./TensorComp/3d/W_Module.jl"); using .ω_module
#include("./Dynamics/Dynamics.jl"); using .Dynamics
include("./Mesh/3d/SLOD_MESH.jl"); using .Mesh3d
include("./Assembly/3d/Assembly.jl"); using .Assemble
include("./Minimization/J_method.jl")
include("Quadrature_3d.jl");
include("./Dynamics/Dynamics.jl");
