"""
    WeakClusterExpansionLattice(max_order, expansion_unit_cell)

An infinite lattice for the weak-embedding cluster expansion up to `max_order`.
"""
struct WeakClusterExpansionLattice <: AbstractClusterExpansionLattice
        max_order::UInt8

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

function WeakClusterExpansionLattice(max_order::Int, expansion_unit_cell::ExpansionUnitCell)
        @assert max_order > 0 "max_order must be a positive integer"
        @assert !is_strong_tiling(expansion_unit_cell) "Strong tiling - use StrongClusterExpansionLattice"

        expansion_coordinates = generate_coordinates(ntuple(_ -> max_order, dimension(expansion_unit_cell)), length(basis_size(expansion_unit_cell)))
        neighbor_list = generate_neighbor_list(expansion_coordinates, expansion_unit_cell)
        center_vertices = ExpansionVertices(find_centers(expansion_coordinates))

        lattice_coordinates = build_lattice_coordinates(max_order, expansion_unit_cell)

        connections_vec, rev_connections, masking_matrix, unique_inds = generate_weak_connections(expansion_coordinates, lattice_coordinates, expansion_unit_cell)

        adj_matrix = generate_adj_matrix_weak(lattice_coordinates, expansion_unit_cell, unique_inds)

        lattice_coordinates = lattice_coordinates[:, unique_inds]
        translation_labels = [expansion_unit_cell.translation_labels[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        site_colors = [expansion_unit_cell.site_colors[col[end-1]][col[end]] for col in eachcol(lattice_coordinates)]
        n_unique = length(unique(Iterators.flatten(expansion_unit_cell.translation_labels)))

        return WeakClusterExpansionLattice(
                UInt8(max_order),
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

