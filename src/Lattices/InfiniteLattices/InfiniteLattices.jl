"""
    AbstractInfiniteLattice <: AbstractLattice

Abstract base type for an infinite lattice used to enumerate clusters.

Subtypes must implement all methods of `AbstractLattice`, plus:
- `n_unique_sites(lattice)` - number of translationally-inequivalent sites in the lattice
"""
abstract type AbstractInfiniteLattice <: AbstractLattice end

"""
    AbstractClusterExpansionLattice <: AbstractInfiniteLattice

Abstract base type for lattices where each expansion vertex represents a cluster of
physical sites (a unit cell), rather than a single site.

Subtypes must implement all methods of `AbstractInfiniteLattice`, plus:
- `n_site_colors(lattice)` - number of distinct site colors in the lattice 
- `connections(lattice)` - mapping from expansion vertices to their constituent lattice vertices
"""
abstract type AbstractClusterExpansionLattice <: AbstractInfiniteLattice end

n_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

centers(lattice::AbstractClusterExpansionLattice) = lattice.centers
max_order(lattice::AbstractClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::AbstractClusterExpansionLattice) = lattice.n_unique
n_site_colors(lattice::AbstractClusterExpansionLattice) = length(unique(Iterators.flatten(lattice.expansion_unit_cell.site_colors)))
neighbors(lattice::AbstractClusterExpansionLattice, vertices::ExpansionVertices) = setdiff(union(ExpansionVertices(), lattice.neighbor_list[vertices]), vertices)
get_coordinates(lattice::AbstractClusterExpansionLattice) = shift_unit_cell(lattice.expansion_unit_cell, lattice.lattice_coordinates)
get_labels(lattice::AbstractClusterExpansionLattice) = lattice.translation_labels
get_site_colors(lattice::AbstractClusterExpansionLattice) = lattice.site_colors
connections(lattice::AbstractClusterExpansionLattice) = lattice.connections
bond_matrix(lattice::AbstractClusterExpansionLattice) = lattice.adj_matrix

include("SiteExpansionLattices.jl")
include("StrongClusterExpansionLattices.jl")
include("WeakClusterExpansionLattices.jl")
