"""
    ExpansionCluster(cluster, lattice)

Cluster from a set of clusters that contains information necessary for summation! and writing to disk.
"""
struct ExpansionCluster <: AbstractExpansionCluster
        vertices::LatticeVertices{Int}
        lattice_constant::Float64
        subgraphs::Vector{UInt}
        weights::Dict{UInt,Float64}
end

function _build_site_expansion_cluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::AbstractLattice)
        subgraphs = Vector{UInt}()
        for subgraph in get_subgraphs(cluster, lattice)
                push!(subgraphs, ghash(clusters, subgraph))
        end

        ExpansionCluster(cluster.vertices, lattice_constant(cluster), subgraphs, Dict{UInt,Float64}(cluster.ghash => 1.0))
end

function _build_expansion_cluster(lattice_vertices::LatticeVertices, cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::AbstractLattice)
        subgraphs = Vector{UInt}()
        sizehint!(subgraphs, length(lattice_vertices) + length(cluster) - 1)

        for lv in lattice_vertices
                push!(subgraphs, ghash(clusters, LatticeVertices(lv)))
        end
        for subgraph in get_subgraphs(cluster, lattice)
                push!(subgraphs, ghash(clusters, subgraph))
        end

        ExpansionCluster(lattice_vertices, lattice_constant(cluster), subgraphs, Dict{UInt,Float64}(cluster.ghash => 1.0))
end

ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::SiteExpansionLattice) =
        _build_site_expansion_cluster(cluster, clusters, lattice)

ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::FiniteLattice) =
        _build_site_expansion_cluster(cluster, clusters, lattice)

ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::Union{StrongClusterExpansionLattice,FiniteStrongClusterExpansionLattice}) =
        _build_expansion_cluster(connections(lattice)[cluster.vertices], cluster, clusters, lattice)

ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::Union{WeakClusterExpansionLattice,FiniteWeakClusterExpansionLattice}) =
        _build_expansion_cluster(just_lattice_vertices(connections(lattice), cluster.vertices), cluster, clusters, lattice)

function ExpansionCluster(lv::Int, single_site_hash::UInt, n_single_site_clusters::Int)
        cluster_lattice_constant = 1 / n_single_site_clusters

        ExpansionCluster(LatticeVertices(lv), cluster_lattice_constant, UInt[], Dict{UInt,Float64}(single_site_hash => 1.0))
end

lattice_constant(cluster::ExpansionCluster) = cluster.lattice_constant
subgraphs(cluster::ExpansionCluster) = cluster.subgraphs
function subtract_subcluster!(cluster::ExpansionCluster, subcluster::ExpansionCluster)
        for (k, v) in subcluster.weights
                cluster.weights[k] = get(cluster.weights, k, 0.0) - v
        end
end
