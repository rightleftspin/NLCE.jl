@testset verbose = true "Clusters" begin

        @testset "Translation clustering (Square)" begin
                lattice = SiteExpansionLattice(4, square_uc)
                clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(clusters, lattice)

                @test length(filter(c -> length(c) == 1, collect(clusters))) == 1
                @test length(filter(c -> length(c) == 2, collect(clusters))) == 2
                @test length(filter(c -> length(c) == 3, collect(clusters))) == 6
                @test length(filter(c -> length(c) == 4, collect(clusters))) == 19
        end

        @testset "Connected clustering (Square)" begin
                lattice = SiteExpansionLattice(4, square_uc)
                clusters = ConnectedClusterSet(lattice)
                clusters_from_lattice!(clusters, lattice)

                @test length(filter(c -> length(c) == 1, collect(clusters))) == 1
                @test length(filter(c -> length(c) == 2, collect(clusters))) == 4
                @test length(filter(c -> length(c) == 3, collect(clusters))) == 18
                @test length(filter(c -> length(c) == 4, collect(clusters))) == 76
        end

        @testset "Symmetric clustering (Square)" begin
                lattice = SiteExpansionLattice(8, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                sym_clusters = SymmetricClusterSet(lattice, :Square)
                clusters_from_clusters!(sym_clusters, trans_clusters)

                @test length(sym_clusters) == 533
        end

        @testset "Isomorphic clustering (Square)" begin
                lattice = SiteExpansionLattice(4, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                order4 = filter(c -> length(c) == 4, collect(iso_clusters))
                @test length(order4) == 3
        end

        @testset "Translation clustering (Pyrochlore Unit Cell)" begin
                lattice = StrongClusterExpansionLattice(3, pyro_exp_uc_uc)
                clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(clusters, lattice)

                @test length(filter(c -> length(c) == 1, collect(clusters))) == 1
                @test length(filter(c -> length(c) == 2, collect(clusters))) == 6
                @test length(filter(c -> length(c) == 3, collect(clusters))) == 50
        end

        @testset "Isomorphic clustering (Pyrochlore Unit Cell)" begin
                lattice = StrongClusterExpansionLattice(4, pyro_exp_uc_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                order3 = filter(c -> length(c) == 3, collect(iso_clusters))
                order4 = filter(c -> length(c) == 4, collect(iso_clusters))
                @test length(order3) == 3
                @test length(order4) == 8
        end

        @testset "Translation clustering (Square Cluster)" begin
                lattice = WeakClusterExpansionLattice(4, square_cluster_uc)
                clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(clusters, lattice)

                @test length(filter(c -> length(c) == 1, collect(clusters))) == 1
                @test length(filter(c -> length(c) == 2, collect(clusters))) == 2
                @test length(filter(c -> length(c) == 3, collect(clusters))) == 6
                @test length(filter(c -> length(c) == 4, collect(clusters))) == 19
        end

        @testset "Symmetric clustering (Square Cluster)" begin
                lattice = WeakClusterExpansionLattice(6, square_cluster_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                sym_clusters = SymmetricClusterSet(lattice, :Square)
                clusters_from_clusters!(sym_clusters, trans_clusters)

                @test length(sym_clusters) == 56
        end

        @testset "Isomorphic clustering (Square Cluster)" begin
                lattice = WeakClusterExpansionLattice(5, square_cluster_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                order3 = filter(c -> length(c) == 3, collect(iso_clusters))
                order4 = filter(c -> length(c) == 4, collect(iso_clusters))
                order5 = filter(c -> length(c) == 5, collect(iso_clusters))
                @test length(order3) == 2
                @test length(order4) == 5
                @test length(order5) == 11
        end
end # Clusters
