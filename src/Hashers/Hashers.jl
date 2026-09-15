"""
    AbstractHasher

Abstract base type for cluster hashing. A hasher maps a vertex set to a graph hash, with two vertex sets receiving the same hash if and only if they are equivalent under the symmetry the hasher preserves.

Subtypes must implement:
- `ghash(hasher, vertices::LatticeVertices)` - hash for a Site Expansions
- `ghash(hasher, vertices::ExpansionVertices)` - hash for a Cluster Expansion 
"""
abstract type AbstractHasher end

ghash(h::AbstractHasher, expansion_vertices::ExpansionVertices) = _NI("ghash")
ghash(h::AbstractHasher, lattice_vertices::LatticeVertices) = _NI("ghash")
n_unique_sites(h::AbstractHasher) = _NI("n_unique_sites")

include("util.jl")
include("TranslationHasher.jl")
include("IsomorphicHasher.jl")
include("SymmetricHasher.jl")
include("ConnectedHasher.jl")

