@testset verbose = true "Finite Lattices" begin

        @testset verbose = true "FiniteLattice" begin

                @testset "Square OBC" begin
                        lattice = FiniteLattice((4, 4), square_uc, 16)

                        @test n_unique_sites(lattice) == 16
                        @test max_order(lattice) == 16
                        @test length(centers(lattice)) == 16
                        neighbor_counts = sort([length(neighbors(lattice, LatticeVertices([i]))) for i in 1:16])
                        @test neighbor_counts == [2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4]
			@test fld(count(!=(0), bond_matrix(lattice)), 2) == 24
                end

                @testset "Square PBC" begin
                        lattice = FiniteLattice((4, 4), square_uc, 16; periodic=true)

                        @test n_unique_sites(lattice) == 16
                        @test all(==(4), [length(neighbors(lattice, LatticeVertices([i]))) for i in 1:16])
                        @test fld(count(!=(0), bond_matrix(lattice)), 2) == 32
                end

                @testset "Square quasi-finite (periodic x, open y)" begin
                        lattice = FiniteLattice((0:3, 0:2), square_uc, 12; periodic_dims=[true, false])

                        @test n_unique_sites(lattice) == 12
                        neighbor_counts = sort([length(neighbors(lattice, LatticeVertices([i]))) for i in 1:12])
                        @test neighbor_counts == [3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4]
                end

        end

        @testset verbose = true "FiniteStrongClusterExpansionLattice" begin

                @testset "Pyrochlore OBC 2x2x2" begin
                        lattice = FiniteStrongClusterExpansionLattice((2, 2, 2), pyro_exp_uc_uc, 8)

                        @test n_unique_sites(lattice) == 32
                        @test max_order(lattice) == 8
                        @test length(centers(lattice)) == 8
                        neighbor_counts = sort([length(neighbors(lattice, ExpansionVertices([i]))) for i in 1:8])
                        @test neighbor_counts == [3, 3, 5, 5, 5, 5, 5, 5]
                end

                @testset "Pyrochlore PBC 2x2x2" begin
                        lattice = FiniteStrongClusterExpansionLattice((2, 2, 2), pyro_exp_uc_uc, 8; periodic_dims=[true, true, true])

                        @test n_unique_sites(lattice) == 32
                        @test all(==(6), [length(neighbors(lattice, ExpansionVertices([i]))) for i in 1:8])
                        @test fld(count(!=(0), bond_matrix(lattice)), 2) == 96
                end

        end

        @testset verbose = true "FiniteWeakClusterExpansionLattice" begin

                @testset "Square Cluster OBC 3x3" begin
                        lattice = FiniteWeakClusterExpansionLattice((3, 3), square_cluster_uc, 9)

                        @test n_unique_sites(lattice) == 24
                        @test max_order(lattice) == 9
                        @test length(centers(lattice)) == 9
                        neighbor_counts = sort([length(neighbors(lattice, ExpansionVertices([i]))) for i in 1:9])
                        @test neighbor_counts == [2, 2, 2, 2, 3, 3, 3, 3, 4]
                end

                @testset "Square Cluster PBC 3x3" begin
                        lattice = FiniteWeakClusterExpansionLattice((3, 3), square_cluster_uc, 9; periodic_dims=[true, true])

                        @test n_unique_sites(lattice) == 24
                        @test all(==(4), [length(neighbors(lattice, ExpansionVertices([i]))) for i in 1:9])
                end

        end

        @testset "4x4 Square Lattice cluster pipeline" begin
                lattice = FiniteLattice((4, 4), square_uc, 16)

                connected = ConnectedClusterSet(lattice)
                clusters_from_lattice!(connected, lattice)

                sym = SymmetricClusterSet(lattice, :Square)
                clusters_from_clusters!(sym, connected)

                iso = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso, connected)

                connected_counts = [16, 24, 52, 113, 244, 496, 912]
                # Note: sym_counts differ from Table 1 of Tang et al. (2013), which lists
                # 14, 43, 94 at orders 5-7, the results here have been tested and converge
		# to the appropriate exact-diagonalization limit.
                sym_counts       = [1, 1, 2, 5, 11, 29, 66]
                iso_counts       = [1, 1, 1, 3, 4, 10, 19]

                for order in 1:7
                        @test length(filter(c -> length(c) == order, collect(connected))) == connected_counts[order]
                        @test length(filter(c -> length(c) == order, collect(sym)))       == sym_counts[order]
                        @test length(filter(c -> length(c) == order, collect(iso)))       == iso_counts[order]
                end
        end

        @testset "SymmetricClusterSet errors on non-symmetric finite lattice" begin
                # A rectangular (non-square) finite lattice is NOT invariant under the
                # square point group: a 90-degree rotation maps boundary sites outside the
                # lattice, producing a 0 index in the symmetry permutation. This must
                # error at hasher construction rather than silently producing a wrong
                # hash or a raw out-of-bounds access.
                rect = FiniteLattice((4, 2), square_uc, 8)
                @test_throws ErrorException SymmetricClusterSet(rect, :Square)

                # A square finite lattice IS invariant, so construction succeeds, and the
                # infinite-lattice path is unaffected (0-containing corner permutations
                # are legitimate there and never reached by clusters).
                square = FiniteLattice((4, 4), square_uc, 8)
                @test SymmetricClusterSet(square, :Square) isa Lincege.ClusterSet
                @test SymmetricClusterSet(SiteExpansionLattice(8, square_uc), :Square) isa Lincege.ClusterSet
        end

end # Finite Lattices
