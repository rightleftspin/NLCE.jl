"""
    Bond(site1, site2, direction, bond_type)

Bond of given type between two sites along the given direction, in primitive lattice vector notation
"""
struct Bond
        site1::Int
        site2::Int
        direction::Vector{Int}
        bond_type::Int
end

neighbor_site(bond::Bond, coordinate::AbstractVector{Int}) = [coordinate[1:end-1] + bond.direction; bond.site2]

function pbc_neighbor_fn(dims::NTuple{N,Int}, periodic_dims::Vector{Bool}=fill(true, N)) where {N}
        (bond, col) -> begin
                spatial = col[1:N] .+ bond.direction
                wrapped = [periodic_dims[i] ? mod(spatial[i], dims[i]) : spatial[i] for i in 1:N]
                [wrapped; col[N+1:end-length(bond.site2)]; bond.site2]
        end
end

function supercell_neighbor_fn(supercell_vecs::Matrix{Int})
        S    = Float64.(supercell_vecs)
        Sinv = inv(S)
        N    = size(supercell_vecs, 1)
        (bond, col) -> begin
                spatial = col[1:N] .+ bond.direction
                wrapped = S * mod.(Sinv * spatial, 1.0)
                [round.(Int, wrapped); col[N+1:end-length(bond.site2)]; bond.site2]
        end
end

"""
    UnitCell(basis, primitive_vectors, bonds, site_colors)

Unit cell containing the given information, generally to be used for a site expansion. Coordinates are written [x1 x2 x3 ...; y1 y2 y3 ...; z1 z2 z3 ...;].
"""
struct UnitCell <: AbstractUnitCell
        basis::Matrix{Float64}
        primitive_vectors::Matrix{Float64}
        site_colors::Vector{Int}
        bonds::Vector{Bond}
end

function UnitCell(basis::AbstractVector{<:AbstractVector{Float64}}, primitive_vectors::AbstractVector{<:AbstractVector{Float64}}, bonds::AbstractVector{Bond}, site_colors::AbstractVector{Int})

        @assert length(site_colors) == length(basis) "Number of site colors must match number of sites in the basis."
        @assert all(length.(primitive_vectors) .== length(primitive_vectors)) "Primitive vectors must form a square matrix"
        @assert all(length.(basis) .== length(basis[1])) "Basis vectors must all have the same dimension"


        UnitCell(hcat(basis...), hcat(primitive_vectors...), site_colors, bonds)
end


basis_size(unit_cell::UnitCell) = size(unit_cell.basis, 2)
dimension(unit_cell::UnitCell) = size(unit_cell.primitive_vectors, 2)
shift_unit_cell(unit_cell::UnitCell, shift_vector::AbstractVector{Int}) = unit_cell.primitive_vectors * shift_vector[1:end-1] + unit_cell.basis[:, shift_vector[end]]
shift_unit_cell(unit_cell::UnitCell, shift_vectors::AbstractMatrix{Int}) = stack(v -> shift_unit_cell(unit_cell, v), eachcol(shift_vectors))
