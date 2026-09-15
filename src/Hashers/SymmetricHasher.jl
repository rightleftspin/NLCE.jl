"""
    SymmetricHasher(lattice, lattice_symmetries)

Hasher that identifies clusters equivalent under the lattice's point-group symmetries.
"""
struct SymmetricHasher{C<:Union{<:AbstractConnections,Nothing},T<:TranslationHasher} <: AbstractHasher
        trans_hasher::T
        permutations::Vector{Vector{Int64}}
        connections::C
end

function SymmetricHasher(lattice::AbstractLattice, lattice_symmetries::Vector{Matrix{Float64}}, connections::Union{<:AbstractConnections,Nothing})
        trans_hasher = TranslationHasher(lattice)
        all_coords = get_coordinates(lattice)
        # `get_permutations` may produce incomplete permutations (entries set to 0)
        # when a symmetry maps a boundary site outside the lattice's coordinate cube.
        # For infinite-lattice use this is safe: clusters are grown from the center
        # and never reach the corners, so no cluster site will land on a zero entry.
        permutations = get_permutations(all_coords, lattice_symmetries)
        _validate_finite_permutations(lattice, permutations)

        SymmetricHasher(
                trans_hasher,
                permutations,
                connections
        )

end

_validate_finite_permutations(::AbstractLattice, ::Vector{Vector{Int64}}) = nothing
function _validate_finite_permutations(::AbstractFiniteLattice, permutations::Vector{Vector{Int64}})
        for perm in permutations
                if any(iszero, perm)
                        error("SymmetricHasher: a lattice point-group symmetry maps a site outside the finite lattice. Use a lattice whose shape is invariant under the given symmetries (e.g. a 4x4 grid for the square lattice).")
                end
        end
        nothing
end

SymmetricHasher(lattice::AbstractLattice, lattice_symmetries::Vector{Matrix{Float64}}) = SymmetricHasher(lattice, lattice_symmetries, nothing)
SymmetricHasher(lattice::AbstractClusterExpansionLattice, lattice_symmetries::Vector{Matrix{Float64}}) = SymmetricHasher(lattice, lattice_symmetries, connections(lattice))
SymmetricHasher(lattice::AbstractFiniteClusterExpansionLattice, lattice_symmetries::Vector{Matrix{Float64}}) = SymmetricHasher(lattice, lattice_symmetries, connections(lattice))

n_unique_sites(h::SymmetricHasher) = n_unique_sites(h.trans_hasher)

function ghash(h::SymmetricHasher, lattice_vertices::LatticeVertices)
        idx = collect(lattice_vertices)
        # minimum over the orbit is a canonical representative
        minimum(ghash_idx(h.trans_hasher, sort(perm[idx])) for perm in h.permutations)
end

ghash(h::SymmetricHasher{StrongClusterConnections}, expansion_vertices::ExpansionVertices) = ghash(h, h.connections[expansion_vertices])

function ghash(h::SymmetricHasher{WeakClusterConnections}, expansion_vertices::ExpansionVertices)
        lattice_vertices, mask = h.connections[expansion_vertices]
        function perm_hash(perm)
                new_lattice_vertices = perm[lattice_vertices]
                sp = sortperm(new_lattice_vertices)
                hm = h.trans_hasher.hashing_matrix[new_lattice_vertices[sp], new_lattice_vertices[sp]]
                hm[mask[sp, sp]] .= 0
                hash(sum(hm, dims=2))
        end
        minimum(perm_hash(perm) for perm in h.permutations)
end
