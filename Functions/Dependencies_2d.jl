using LinearAlgebra, SparseArrays, IterativeSolvers,JLD, MKLSparse


include("./TensorComp/2d/W_Module.jl"); using .ω_module
#include("./Dynamics/Dynamics.jl"); using .Dynamics
include("./Mesh/2d/SLOD_MESH.jl"); using .SLOD_MESH
include("./Assembly/2d/Assembly.jl"); using .Assemble
include("./Minimization/J_method.jl")
include("./Minimization/Sobolev_Gradient.jl")
include("Quadrature_2d.jl");
include("./Dynamics/Dynamics.jl");
