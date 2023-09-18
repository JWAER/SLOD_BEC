module Mesh3d

export Mesh

include("3dmesh.jl");

struct Mesh
    H::Float64
    p::Matrix{Float64}
    t::Matrix{Int}
    box_size::Float64
    CoarseNodes::Vector{Int}
    CoarseNodalLayer::Matrix{Int}
    Node2Tetrahedra::Matrix{Int}
    Faces::Vector{Int}
    Faces_dofs::Matrix{Int}
    bdry::Vector{Int}
    dofs::Vector{Int}
end



struct Mesh_refined
    H::Float64
    p::Matrix{Float64}
    t::Matrix{Int}
    t_H ::Matrix{Int}
    box_size::Float64
    CoarseNodes::Vector{Int}
    CoarseNodalLayer::Matrix{Int}
    CoarseNode2CoarseSimplices::Matrix{Int}
    FineSimplicesInCoarse::Vector{UnitRange{Int64}}
    Faces::Vector{Int}
    Faces_dofs::Matrix{Int}
    bdry::Vector{Int}
    dofs::Vector{Int}
end



end
