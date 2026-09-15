function add_cluster!(cluster_set::AbstractClusterSet{C,H}, c::AbstractCluster, graph_hash::UInt) where {C,H}
        old = get(cluster_set, graph_hash, nothing)
        push!(cluster_set, C(c.vertices, isnothing(old) ? c.lattice_constant : c.lattice_constant + old.lattice_constant, graph_hash))
end

add_cluster!(cluster_set::AbstractClusterSet, c::AbstractCluster) = add_cluster!(cluster_set, c, ghash(cluster_set, c))

"""
    clusters_from_clusters!(new_clusters, old_clusters)

Populates new_clusters from the clusters inside old_clusters using the hasher inside the new_clusters ClusterSet
"""
function clusters_from_clusters!(new_clusters::AbstractClusterSet{C,H}, old_clusters::AbstractClusterSet) where {C<:AbstractCluster,H}
        old_vec = collect(old_clusters)
        hashes = Vector{UInt}(undef, length(old_vec))
        @threads for i in eachindex(old_vec)
                hashes[i] = ghash(new_clusters, old_vec[i])
        end
        for i in eachindex(old_vec)
                add_cluster!(new_clusters, old_vec[i], hashes[i])
        end
        new_clusters
end

"""
    clusters_from_lattice!(clusters, lattice; spawn_depth=3)

Generates all clusters from an infinite lattice up till the given max_order, populates clusters with all the corresponding clusters reduced by the hashing function
"""
function clusters_from_lattice!(clusters::AbstractClusterSet{C,H}, lattice::AbstractLattice; spawn_depth::Int=3) where {C<:AbstractCluster,H}
        max_depth = max_order(lattice)
        center_vertices = centers(lattice)
        roots = [Cluster(typeof(center_vertices)(center), clusters, lattice) for center in center_vertices]
        vlock = ReentrantLock()

        function try_mark(cluster::C)
                already = lock(vlock) do
                        already_present = cluster in clusters
                        if !already_present
                                push!(clusters, cluster)
                        end
                        already_present
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster::C, depth::Int)
                if !try_mark(cluster)
                        return
                end

                if depth < spawn_depth
                        for v in neighbors(lattice, cluster.vertices)
                                dfs(Cluster(union(cluster.vertices, typeof(center_vertices)(v)), clusters, lattice), depth + 1)
                        end
                else
                        @sync begin
                                first = true
                                for v in neighbors(lattice, cluster.vertices)
                                        neighbor_cluster = Cluster(union(cluster.vertices, typeof(center_vertices)(v)), clusters, lattice)

                                        if first
                                                # Run the first neighbor inline to avoid unnecessary task allocation,
                                                # then spawn the remaining neighbors concurrently.
                                                dfs(neighbor_cluster, depth + 1)
                                                first = false
                                        else
                                                @spawn dfs(neighbor_cluster, depth + 1)
                                        end
                                end
                        end
                end
        end

        for root in roots
                dfs(root, 0)
        end

        clusters
end
