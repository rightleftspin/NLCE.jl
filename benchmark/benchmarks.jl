using Lincege
using BenchmarkTools

SUITE = BenchmarkGroup()

square_uc = Lincege.square_unit_cell
lattice_sq = SiteExpansionLattice(10, square_uc)
trans_sq = TranslationClusterSet(lattice_sq)
clusters_from_lattice!(trans_sq, lattice_sq)
iso_sq = IsomorphicClusterSet(lattice_sq)
clusters_from_clusters!(iso_sq, trans_sq)

SUITE["square"] = BenchmarkGroup()
SUITE["square"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_sq)
        clusters_from_lattice!(cs, $lattice_sq)
end
SUITE["square"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_sq)
        clusters_from_clusters!(iso, $trans_sq)
end
SUITE["square"]["clusters_from_clusters_symmetric"] = @benchmarkable begin
        sym = SymmetricClusterSet($lattice_sq, :Square)
        clusters_from_clusters!(sym, $trans_sq)
end
SUITE["square"]["Expansion"] = @benchmarkable Expansion($iso_sq, $lattice_sq)
SUITE["square"]["summation"] = @benchmarkable begin
        e = Expansion($iso_sq, $lattice_sq)
        summation!(e, 10)
end

kagome_uc = Lincege.kagome_unit_cell
lattice_kag = SiteExpansionLattice(4, kagome_uc)
trans_kag = TranslationClusterSet(lattice_kag)
clusters_from_lattice!(trans_kag, lattice_kag)
iso_kag = IsomorphicClusterSet(lattice_kag)
clusters_from_clusters!(iso_kag, trans_kag)

SUITE["kagome"] = BenchmarkGroup()
SUITE["kagome"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_kag)
        clusters_from_lattice!(cs, $lattice_kag)
end
SUITE["kagome"]["summation"] = @benchmarkable begin
        e = Expansion($iso_kag, $lattice_kag)
        summation!(e, 4)
end

# Pyrochlore Unit Cell - StrongClusterExpansionLattice
pyro_exp_uc_uc = Lincege.pyrochlore_expansion_unit_cell
lattice_pyro = StrongClusterExpansionLattice(5, pyro_exp_uc_uc)
trans_pyro = TranslationClusterSet(lattice_pyro)
clusters_from_lattice!(trans_pyro, lattice_pyro)
iso_pyro = IsomorphicClusterSet(lattice_pyro)
clusters_from_clusters!(iso_pyro, trans_pyro)

SUITE["pyrochlore"] = BenchmarkGroup()
SUITE["pyrochlore"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_pyro)
        clusters_from_lattice!(cs, $lattice_pyro)
end
SUITE["pyrochlore"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_pyro)
        clusters_from_clusters!(iso, $trans_pyro)
end
SUITE["pyrochlore"]["summation"] = @benchmarkable begin
        e = Expansion($iso_pyro, $lattice_pyro)
        summation!(e, 5)
end

# Square Cluster - WeakClusterExpansionLattice
square_cluster_uc = Lincege.square_cluster_expansion_unit_cell
lattice_sq_cluster = WeakClusterExpansionLattice(4, square_cluster_uc)
trans_sq_cluster = TranslationClusterSet(lattice_sq_cluster)
clusters_from_lattice!(trans_sq_cluster, lattice_sq_cluster)
iso_sq_cluster = IsomorphicClusterSet(lattice_sq_cluster)
clusters_from_clusters!(iso_sq_cluster, trans_sq_cluster)

SUITE["square_cluster"] = BenchmarkGroup()
SUITE["square_cluster"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_sq_cluster)
        clusters_from_lattice!(cs, $lattice_sq_cluster)
end
SUITE["square_cluster"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_sq_cluster)
        clusters_from_clusters!(iso, $trans_sq_cluster)
end
SUITE["square_cluster"]["clusters_from_clusters_symmetric"] = @benchmarkable begin
        sym = SymmetricClusterSet($lattice_sq_cluster, :Square)
        clusters_from_clusters!(sym, $trans_sq_cluster)
end
SUITE["square_cluster"]["summation"] = @benchmarkable begin
        e = Expansion($iso_sq_cluster, $lattice_sq_cluster)
        summation!(e, 4)
end

# Finite lattices
SUITE["finite"] = BenchmarkGroup()

# 4x4 open-boundary square lattice (site expansion)
SUITE["finite"]["square_connected"] = @benchmarkable begin
        lat = FiniteLattice((4, 4), $square_uc, 16)
        cs = ConnectedClusterSet(lat)
        clusters_from_lattice!(cs, lat)
end

# 2x2x2 open-boundary pyrochlore (3D strong cluster expansion)
SUITE["finite"]["pyrochlore_strong"] = @benchmarkable begin
        lat = FiniteStrongClusterExpansionLattice((2, 2, 2), $pyro_exp_uc_uc, 8)
        cs = ConnectedClusterSet(lat)
        clusters_from_lattice!(cs, lat)
end
