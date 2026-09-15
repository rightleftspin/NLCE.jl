"""
    FiniteStrongClusterExpansionLattice(dims, expansion_unit_cell, max_order; periodic_dims=fill(false, N))

A finite lattice for a *strong-embedding* cluster expansion of size `dims` unit cells. Each
expansion vertex is a cluster (unit cell) of physical sites.

Pass `periodic_dims` to choose which dimensions use periodic boundary conditions; by default all
boundaries are open.
"""
struct FiniteStrongClusterExpansionLattice <: AbstractFiniteClusterExpansionLattice
        max_order::Int

        expansion_unit_cell::ExpansionUnitCell
        centers::ExpansionVertices
        neighbor_list::Vector{ExpansionVertices{Int}}

        lattice_coordinates::Matrix{Int}
        translation_labels::Vector{Int}
        site_colors::Vector{Int}
        adj_matrix::Matrix{Int}

        connections::StrongClusterConnections
        n_unique::Int
end

function FiniteStrongClusterExpansionLattice(dims::NTuple{N,Int}, expansion_unit_cell::ExpansionUnitCell, max_order::Int; periodic_dims::Vector{Bool}=fill(false, N)) where {N}
        @assert is_strong_tiling(expansion_unit_cell) "Not a strong tiling - use FiniteWeakClusterExpansionLattice"
        expansion_coordinates = generate_finite_coordinates(dims, length(basis_size(expansion_unit_cell)))
        center_vertices = ExpansionVertices(1:size(expansion_coordinates, 2))

        lattice_coordinates = build_lattice_coordinates_finite(dims, expansion_unit_cell)
        translation_labels = [expansion_unit_cell.translation_labels[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        site_colors = [expansion_unit_cell.site_colors[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        connections_vec = generate_strong_connections(expansion_coordinates, lattice_coordinates)

        if any(periodic_dims)
                nfn           = pbc_neighbor_fn(dims, periodic_dims)
                neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell; neighbor_fn=nfn)
                adj_matrix    = generate_adj_matrix(lattice_coordinates, expansion_unit_cell; neighbor_fn=nfn)
        else
                neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell)
                adj_matrix    = generate_adj_matrix(lattice_coordinates, expansion_unit_cell)
        end

        n_unique = prod(dims) * sum(basis_size(expansion_unit_cell))
        return FiniteStrongClusterExpansionLattice(
                max_order,
                expansion_unit_cell,
                center_vertices,
                neighbor_list,
                lattice_coordinates,
                translation_labels,
                site_colors,
                adj_matrix,
                StrongClusterConnections(connections_vec),
                n_unique
        )
end

"""
    FiniteWeakClusterExpansionLattice(dims, expansion_unit_cell, max_order; periodic_dims=fill(false, N))

A finite lattice for a *weak-embedding* cluster expansion of size `dims` unit cells. Each
expansion vertex is a cluster of physical sites.

Pass `periodic_dims` to choose which dimensions use periodic boundary conditions; by default all
boundaries are open.
"""
struct FiniteWeakClusterExpansionLattice <: AbstractFiniteClusterExpansionLattice
        max_order::Int

        expansion_unit_cell::ExpansionUnitCell
        centers::ExpansionVertices
        neighbor_list::Vector{ExpansionVertices{Int}}

        lattice_coordinates::Matrix{Int}
        translation_labels::Vector{Int}
        site_colors::Vector{Int}
        adj_matrix::Matrix{Int}

        connections::WeakClusterConnections
        n_unique::Int
end

function FiniteWeakClusterExpansionLattice(dims::NTuple{N,Int}, expansion_unit_cell::ExpansionUnitCell, max_order::Int; periodic_dims::Vector{Bool}=fill(false, N)) where {N}
        @assert !is_strong_tiling(expansion_unit_cell) "Strong tiling - use FiniteStrongClusterExpansionLattice"
        expansion_coordinates = generate_finite_coordinates(dims, length(basis_size(expansion_unit_cell)))
        center_vertices = ExpansionVertices(1:size(expansion_coordinates, 2))

        lattice_coordinates = build_lattice_coordinates_finite(dims, expansion_unit_cell)

        connections_vec, rev_connections, masking_matrix, unique_inds = generate_weak_connections(expansion_coordinates, lattice_coordinates, expansion_unit_cell)

        if any(periodic_dims)
                nfn           = pbc_neighbor_fn(dims, periodic_dims)
                neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell; neighbor_fn=nfn)
                adj_matrix    = generate_adj_matrix_weak(lattice_coordinates, expansion_unit_cell, unique_inds; neighbor_fn=nfn)
        else
                neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell)
                adj_matrix    = generate_adj_matrix_weak(lattice_coordinates, expansion_unit_cell, unique_inds)
        end

        lattice_coordinates = lattice_coordinates[:, unique_inds]
        translation_labels = [expansion_unit_cell.translation_labels[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        site_colors = [expansion_unit_cell.site_colors[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        n_unique = length(unique_inds)

        return FiniteWeakClusterExpansionLattice(
                max_order,
                expansion_unit_cell,
                center_vertices,
                neighbor_list,
                lattice_coordinates,
                translation_labels,
                site_colors,
                adj_matrix,
                WeakClusterConnections(connections_vec, rev_connections, masking_matrix, basis_size(expansion_unit_cell)),
                n_unique
        )
end
