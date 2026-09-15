struct ClusterSet{C<:AbstractCluster,H<:AbstractHasher} <: AbstractClusterSet{C,H}
        clusters::Dict{UInt,C}
        hasher::H
end

"""
    TranslationClusterSet(lattice)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves translational invariance
"""
function TranslationClusterSet(lattice::AbstractLattice)
        C = Cluster{typeof(centers(lattice))}
        ClusterSet{C,TranslationHasher}(
                Dict{UInt,C}(),
                TranslationHasher(lattice)
        )
end

"""
    IsomorphicClusterSet(lattice)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves invariance under graph isomorphisms
"""
function IsomorphicClusterSet(lattice::AbstractLattice)
        C = Cluster{typeof(centers(lattice))}
        ClusterSet{C,IsomorphicHasher}(
                Dict{UInt,C}(),
                IsomorphicHasher(lattice)
        )
end

"""
    SymmetricClusterSet(lattice, symmetries)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves invariance under the given point group symmetries of the lattice.
"""
function SymmetricClusterSet(lattice::AbstractLattice, symmetries::Vector{Matrix{Float64}})
        C = Cluster{typeof(centers(lattice))}
        ClusterSet{C,SymmetricHasher}(
                Dict{UInt,C}(),
                SymmetricHasher(lattice, symmetries)
        )
end

SymmetricClusterSet(lattice::AbstractLattice, lattice_type::Symbol) = SymmetricClusterSet(lattice, all_lattice_symmetries[lattice_type])

"""
    ConnectedClusterSet(lattice)

ClusterSet constructor that initializes a ClusterSet with a hasher that treats each distinct set of vertices as unique.
"""
function ConnectedClusterSet(lattice::AbstractLattice)
        C = Cluster{typeof(centers(lattice))}
        ClusterSet{C,ConnectedHasher}(
                Dict{UInt,C}(),
                ConnectedHasher(lattice)
        )
end

Base.length(cluster_set::ClusterSet) = length(cluster_set.clusters)
Base.in(c::C, cluster_set::ClusterSet{C,H}) where {C<:AbstractCluster,H} = haskey(cluster_set.clusters, c.ghash)
Base.iterate(cluster_set::ClusterSet) = iterate(values(cluster_set.clusters))
Base.iterate(cluster_set::ClusterSet, state) = iterate(values(cluster_set.clusters), state)
Base.push!(cluster_set::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = (cluster_set.clusters[c.ghash] = c)
Base.pop!(cluster_set::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = pop!(cluster_set.clusters, c.ghash)
Base.get(cluster_set::ClusterSet{C,H}, graph_hash::UInt, default) where {C<:AbstractCluster,H} = get(cluster_set.clusters, graph_hash, default)
Base.sort(cluster_set::ClusterSet) = sort(collect(values(cluster_set.clusters)), by=length)

ghash(cluster_set::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = ghash(cluster_set.hasher, c.vertices)
ghash(cluster_set::ClusterSet, vertices::AbstractVertices) = ghash(cluster_set.hasher, vertices)
n_unique_sites(cluster_set::ClusterSet) = n_unique_sites(cluster_set.hasher)

Base.show(io::IO, cluster_set::ClusterSet{C,H}) where {C,H} =
        print(io, "ClusterSet{$(nameof(H))} with $(length(cluster_set)) clusters")
