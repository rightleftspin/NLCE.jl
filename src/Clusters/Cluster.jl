struct Cluster{V<:AbstractVertices} <: AbstractCluster
        vertices::V
        lattice_constant::Float64
        ghash::UInt64
end

function Cluster(vertices::AbstractVertices, cluster_set::AbstractClusterSet, lattice::AbstractLattice)
        Cluster(vertices, 1 / n_unique_sites(lattice), ghash(cluster_set, vertices))
end

Base.length(c::Cluster) = length(c.vertices)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
lattice_constant(c::Cluster) = c.lattice_constant
