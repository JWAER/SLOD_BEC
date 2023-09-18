module ω_module


using LinearAlgebra

export Tensor 


	
#------------------------------------------------------------------------------
struct Tensor
	Iptr ::Vector{Int}
	J ::Vector{Int}
	K ::Vector{Int}
	Val ::Vector{Float64}
end
#------------------------------------------------------------------------------

include("Compute_W.jl");
include("Compute_Wh.jl");
include("../PreAllocateW.jl");
include("../../Assembly/3d/Basis.jl");


end
