"""
    TranslationHasher(lattice)

Hasher that identifies clusters equivalent under lattice translations.
"""
struct TranslationHasher{C<:Union{<:AbstractConnections,Nothing}} <: AbstractHasher
        hashing_matrix::Matrix{Int}
        connections::C
        n_unique_sites::Int
end

function TranslationHasher(lattice::AbstractLattice, connections::Union{<:AbstractConnections,Nothing})
        pwd_matrix = pairwise_direction(get_coordinates(lattice))
        hashing_matrix, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
        diag_matrix = diagm(get_labels(lattice) .+ max_dir)
        TranslationHasher(2 .^ (hashing_matrix + diag_matrix), connections, n_unique_sites(lattice))
end

TranslationHasher(lattice::AbstractLattice) = TranslationHasher(lattice, nothing)
TranslationHasher(lattice::AbstractClusterExpansionLattice) = TranslationHasher(lattice, connections(lattice))
TranslationHasher(lattice::AbstractFiniteClusterExpansionLattice) = TranslationHasher(lattice, connections(lattice))

n_unique_sites(h::TranslationHasher) = h.n_unique_sites

ghash_idx(h::TranslationHasher, idx::AbstractVector{Int}) =
        hash(sum(@view(h.hashing_matrix[idx, idx]), dims=2))
ghash(h::TranslationHasher, lattice_vertices::LatticeVertices) = ghash_idx(h, collect(lattice_vertices))
ghash(h::TranslationHasher{StrongClusterConnections}, expansion_vertices::ExpansionVertices) = ghash(h, union(LatticeVertices(), h.connections[expansion_vertices]))

function ghash(h::TranslationHasher{WeakClusterConnections}, expansion_vertices::ExpansionVertices)
        lattice_vertices, mask = h.connections[expansion_vertices]
        hm = h.hashing_matrix[lattice_vertices, lattice_vertices]
        hm[mask] .= 0
        hash(sum(hm, dims=2))
end
