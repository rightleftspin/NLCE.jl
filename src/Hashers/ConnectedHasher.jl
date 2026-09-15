"""
    ConnectedHasher(lattice)

Hasher that identifies clusters as equivalent if they are the exact same graphs.
"""
struct ConnectedHasher{C<:Union{<:AbstractConnections,Nothing}} <: AbstractHasher
        connections::C
        n_unique_sites::Int
end

function ConnectedHasher(lattice::AbstractLattice, connections::Union{<:AbstractConnections,Nothing})
        ConnectedHasher(connections, n_unique_sites(lattice))
end

ConnectedHasher(lattice::AbstractLattice) = ConnectedHasher(lattice, nothing)
ConnectedHasher(lattice::AbstractClusterExpansionLattice) = ConnectedHasher(lattice, connections(lattice))
ConnectedHasher(lattice::AbstractFiniteClusterExpansionLattice) = ConnectedHasher(lattice, connections(lattice))

n_unique_sites(h::ConnectedHasher) = h.n_unique_sites
        
ghash(h::ConnectedHasher, lattice_vertices::LatticeVertices) = hash(lattice_vertices)
ghash(h::ConnectedHasher{StrongClusterConnections}, expansion_vertices::ExpansionVertices) = hash(union(LatticeVertices(), h.connections[expansion_vertices]))
ghash(h::ConnectedHasher{WeakClusterConnections}, expansion_vertices::ExpansionVertices) = hash(just_lattice_vertices(h.connections, expansion_vertices))
