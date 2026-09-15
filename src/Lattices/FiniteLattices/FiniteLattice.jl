"""
    FiniteLattice(dims, unit_cell, max_order; periodic=false)
    FiniteLattice(dim_specs, unit_cell, max_order; periodic_dims=fill(false, N))
    FiniteLattice(uc_positions, unit_cell, max_order; supercell_vecs=nothing)

A finite lattice built from `unit_cell` for a site expansion.

Boundary conditions:
- `periodic=true` gives full periodic boundary conditions based on the primitive lattice vectors
- The `dim_specs` 
- The `uc_positions` constructor builds a lattice from an explicit list of unit-cell positions, 
  `supercell_vecs` is then needed to be able to determine periodicity along the given axes in
  primitive lattice coordinates.
"""
struct FiniteLattice <: AbstractFiniteLattice
        max_order::Int
        unit_cell::UnitCell
        coordinates::Matrix{Int}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{LatticeVertices{Int}}
end

function FiniteLattice(dims::NTuple{N,Int}, unit_cell::UnitCell, max_order::Int; periodic::Bool=false) where {N}
        coordinates = generate_finite_coordinates(dims, basis_size(unit_cell))
        nfn = periodic ? pbc_neighbor_fn(dims) : neighbor_site
        adj_matrix    = generate_adj_matrix(coordinates, unit_cell; neighbor_fn=nfn)
        neighbor_list = generate_neighbor_list(coordinates, unit_cell; neighbor_fn=nfn)
        FiniteLattice(max_order, unit_cell, coordinates, adj_matrix, neighbor_list)
end

function FiniteLattice(uc_positions::Vector{Vector{Int}}, unit_cell::UnitCell, max_order::Int;
                       supercell_vecs::Union{Vector{Vector{Int}},Nothing}=nothing)
        n = basis_size(unit_cell)
        coordinates = hcat([vcat(pos, i) for pos in uc_positions for i in 1:n]...)
        nfn = isnothing(supercell_vecs) ? neighbor_site : supercell_neighbor_fn(hcat(supercell_vecs...))
        adj_matrix    = generate_adj_matrix(coordinates, unit_cell; neighbor_fn=nfn)
        neighbor_list = generate_neighbor_list(coordinates, unit_cell; neighbor_fn=nfn)
        FiniteLattice(max_order, unit_cell, coordinates, adj_matrix, neighbor_list)
end

function FiniteLattice(dim_specs::NTuple{N,Union{Int,UnitRange{Int}}}, unit_cell::UnitCell, max_order::Int; periodic_dims::Vector{Bool}=fill(false, N)) where {N}
        coordinates = generate_coordinates(dim_specs, basis_size(unit_cell))
        dims        = ntuple(i -> length(dim_specs[i] isa Int ? (-dim_specs[i]:dim_specs[i]) : dim_specs[i]), N)
        nfn = any(periodic_dims) ? pbc_neighbor_fn(dims, periodic_dims) : neighbor_site
        adj_matrix    = generate_adj_matrix(coordinates, unit_cell; neighbor_fn=nfn)
        neighbor_list = generate_neighbor_list(coordinates, unit_cell; neighbor_fn=nfn)
        FiniteLattice(max_order, unit_cell, coordinates, adj_matrix, neighbor_list)
end

centers(lattice::FiniteLattice)        = LatticeVertices(1:size(lattice.coordinates, 2))
max_order(lattice::FiniteLattice)      = lattice.max_order
n_unique_sites(lattice::FiniteLattice) = size(lattice.coordinates, 2)

neighbors(lattice::FiniteLattice, vertices::LatticeVertices) =
        setdiff(union(LatticeVertices(), lattice.neighbor_list[vertices]), vertices)

get_coordinates(lattice::FiniteLattice)  = shift_unit_cell(lattice.unit_cell, lattice.coordinates)
get_labels(lattice::FiniteLattice)       = lattice.coordinates[end, :]
get_site_colors(lattice::FiniteLattice)  = lattice.unit_cell.site_colors[lattice.coordinates[end, :]]
bond_matrix(lattice::FiniteLattice)      = lattice.adj_matrix
