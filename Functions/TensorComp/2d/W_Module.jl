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
include("../PreAllocateW.jl");
include("../../Assembly/2d/Basis.jl");


end
