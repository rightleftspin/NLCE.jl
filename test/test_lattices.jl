@testset verbose = true "Lattices" begin

        @testset "Bond Type Correctness" begin
                sq_lat = SiteExpansionLattice(2, square_uc)
                @test sort(unique(filter(!=(0), bond_matrix(sq_lat)))) == [1]

                ss_lat = SiteExpansionLattice(2, ss_uc)
                @test sort(unique(filter(!=(0), bond_matrix(ss_lat)))) == [1, 2, 3]
        end

        @testset "is_strong_tiling" begin
                @test_throws AssertionError StrongClusterExpansionLattice(2, square_cluster_uc)
                @test_throws AssertionError WeakClusterExpansionLattice(2, pyro_exp_uc_uc)
        end

end # Lattices
