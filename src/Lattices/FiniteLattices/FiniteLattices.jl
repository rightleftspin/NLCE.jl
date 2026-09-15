"""
    AbstractFiniteLattice <: AbstractLattice

Abstract base type for a finite lattice used to enumerate clusters.

Subtypes must implement all methods of `AbstractLattice`, plus:
- `n_unique_sites(lattice)` - total number of physical sites in the lattice
"""
abstract type AbstractFiniteLattice <: AbstractLattice end

"""
    AbstractFiniteClusterExpansionLattice <: AbstractFiniteLattice

Abstract base type for finite lattices where each expansion vertex represents a cluster of
physical sites rather than a single site.

Subtypes must implement all methods of `AbstractFiniteLattice`, plus:
- `n_site_colors(lattice)` - number of distinct site colors in the lattice
- `connections(lattice)` - mapping from expansion vertices to their constituent lattice vertices
"""
abstract type AbstractFiniteClusterExpansionLattice <: AbstractFiniteLattice end

n_unique_sites(lattice::AbstractFiniteLattice) = _NI("n_unique_sites")

centers(lattice::AbstractFiniteClusterExpansionLattice) = lattice.centers
max_order(lattice::AbstractFiniteClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::AbstractFiniteClusterExpansionLattice) = lattice.n_unique
n_site_colors(lattice::AbstractFiniteClusterExpansionLattice) = length(unique(Iterators.flatten(lattice.expansion_unit_cell.site_colors)))
neighbors(lattice::AbstractFiniteClusterExpansionLattice, vertices::ExpansionVertices) = setdiff(union(ExpansionVertices(), lattice.neighbor_list[vertices]), vertices)
get_coordinates(lattice::AbstractFiniteClusterExpansionLattice) = shift_unit_cell(lattice.expansion_unit_cell, lattice.lattice_coordinates)
get_labels(lattice::AbstractFiniteClusterExpansionLattice) = lattice.translation_labels
get_site_colors(lattice::AbstractFiniteClusterExpansionLattice) = lattice.site_colors
connections(lattice::AbstractFiniteClusterExpansionLattice) = lattice.connections
bond_matrix(lattice::AbstractFiniteClusterExpansionLattice) = lattice.adj_matrix

include("FiniteLattice.jl")
include("FiniteClusterExpansionLattices.jl")
