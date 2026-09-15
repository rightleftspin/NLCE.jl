"""
    AbstractCluster

Abstract base type for a single cluster - a connected subgraph of the lattice with an
associated lattice coefficient used in the NLCE summation.

Subtypes must implement:
- `Base.length(c)` - number of sites (vertices) in the cluster
- `Base.hash(c, h)` - the cluster's graph hash (`ghash`)
- `lattice_constant` - the cluster's lattice constant
"""
abstract type AbstractCluster end

"""
    AbstractClusterSet{C<:AbstractCluster, H<:AbstractHasher}

Abstract base type for a collection of unique clusters sharing a common hasher.
Clusters that are equivalent under the hasher's symmetry are merged

Subtypes must implement:
- `Base.length(cluster_set)` - number of stored clusters
- `Base.in(c, cluster_set)` - membership test
- `Base.iterate(cluster_set)` / `Base.iterate(cluster_set, state)` - iteration over clusters
- `Base.push!(cluster_set, c)` - add a cluster
- `Base.pop!(cluster_set, c)` - remove and return a cluster
- `Base.get(cluster_set, ghash, default)` - look up a stored cluster by its graph hash
- `ghash(cluster_set, c)` / `ghash(cluster_set, vertices)` - delegate to the hasher
"""
abstract type AbstractClusterSet{C<:AbstractCluster,H<:AbstractHasher} end

# Cluster Methods
Base.length(c::AbstractCluster) = _NI("Base.length")
Base.hash(c::AbstractCluster, h::UInt) = _NI("Base.hash")
lattice_constant(c::AbstractCluster) = _NI("lattice_constant")

Base.isequal(c1::C, c2::C) where {C<:AbstractCluster} = c1 == c2
Base.:(==)(c1::C, c2::C) where {C<:AbstractCluster} = (hash(c1) == hash(c2))

# Cluster Set Methods
Base.length(cluster_set::AbstractClusterSet)::Int = _NI("Base.length")
Base.in(cluster::C, cluster_set::AbstractClusterSet{C,H}) where {C<:AbstractCluster,H} = _NI("Base.in")
Base.iterate(cluster_set::AbstractClusterSet) = _NI("Base.iterate")
Base.iterate(cluster_set::AbstractClusterSet, state) = _NI("Base.iterate")
Base.push!(cluster_set::AbstractClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = _NI("Base.push!")
Base.pop!(cluster_set::AbstractClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = _NI("Base.pop!")
Base.get(cluster_set::AbstractClusterSet, graph_hash::UInt, default) = _NI("Base.get")
ghash(cluster_set::AbstractClusterSet, c::AbstractCluster) = _NI("ghash")
ghash(cluster_set::AbstractClusterSet, vertices::AbstractVertices) = _NI("ghash")
n_unique_sites(cluster_set::AbstractClusterSet) = _NI("n_unique_sites")

include("util.jl")
include("Cluster.jl")
include("ClusterSets.jl")
