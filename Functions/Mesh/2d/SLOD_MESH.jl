#------------------------------------------------------------------------------
module SLOD_MESH

export Mesh


include("GenerateMesh2D_P3.jl");
struct Mesh
    H::Float64
    t_H::Matrix{Int}
    Nh::Int
    N :: Int
    box :: Int
    t_h::Matrix{Int}
    p::Matrix{Float64}
    loc2glob::Matrix{Int}
    CoarseNodes::Vector{Int} #1:(N+1)^2 --> (1:(3N*Nh+1)^2)
    CoarseNode2CoarseSimplices::Matrix{Int}
    FineNode2FineSimplices::Matrix{Int}
    NodalLayer::Matrix{Int} #contains (coarse) nodal layer for each node, ell-dependent 
    FineSimplicesInCoarse::Vector{UnitRange{Int64}} # each coarse simplex i is subdivided into fine simplices FineSimplices[i]
    bdry::Vector{Int} #contains boundary nodes 
end
#------------------------------------------------------------------------------



end
