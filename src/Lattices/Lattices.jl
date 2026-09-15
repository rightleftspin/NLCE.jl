"""
    AbstractLattice

Abstract base type for a lattice geometry.

Subtypes must implement:
- `centers(lattice)` - vertex set(s) used as DFS roots when generating clusters
- `max_order(lattice)` - maximum cluster order (number of sites) to generate
- `neighbors(lattice, vertices)` - vertex indices adjacent to the vertex set `vertices` but not in it
- `get_coordinates(lattice)` - Cartesian coordinate matrix (dim x n_sites)
- `get_labels(lattice)` - translation labels for each site
- `get_site_colors(lattice)` - site colors for each site
- `bond_matrix(lattice)` - integer adjacency matrix encoding bond types
"""
abstract type AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vertices::AbstractVertices) = _NI("neighbors")
get_coordinates(lattice::AbstractLattice) = _NI("get_coordinates")
get_labels(lattice::AbstractLattice) = _NI("get_labels")
get_site_colors(lattice::AbstractLattice) = _NI("get_site_colors")
bond_matrix(lattice::AbstractLattice) = _NI("bond_matrix")

include("Connections.jl")
include("util.jl")
include("InfiniteLattices/InfiniteLattices.jl")
include("FiniteLattices/FiniteLattices.jl")
