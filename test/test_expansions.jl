@testset verbose = true "Expansions" begin

        @testset "Square Lattice pipeline" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set(values(weights(expansion, 1))) == Set((1.0))
                @test Set(values(weights(expansion, 2))) == Set((-4.0, 2.0))
                @test Set(values(weights(expansion, 3))) == Set((6.0, -12.0, 6.0))

        end

        @testset "Kagome Lattice pipeline" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, kagome_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((-4.0, 2.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((6.0, -10.0, round(2 / 3, digits=6), 4.0))
        end

        @testset "Pyrochlore Lattice Unit Cell pipeline" begin
                m_order = 3
                lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((0.25, -1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((-3.0, 1.5, 0.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 4))]) == Set((16.5, -27.0, 10.5, 1.0, 1.0, 0.0))
        end

        @testset "Square Lattice Cluster pipeline" begin
                m_order = 3
                lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((0.5, -2.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((-2.0, 1.0, 1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 4))]) == Set((3.0, -6.0, 2.0, 1.0, 0.0))
        end

        @testset "write_to_json" begin
                expected_keys = ["bonds", "cluster_hash", "coordinates",
                        "n_sites", "order", "site_colors", "weights"]

                @testset "SiteExpansionLattice" begin
                        m_order = 3
                        lattice = SiteExpansionLattice(m_order, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["n_sites"] == 1, data)
                                @test length(order1) == 1
                                @test length(order1[1]["bonds"]) == 0
                                @test order1[1]["weights"] == [1, -4, 6]

                                order2 = filter(e -> e["n_sites"] == 2, data)
                                @test length(order2) == 1
                                @test length(order2[1]["bonds"]) == 1
                                @test order2[1]["bonds"][1][end] == 1
                                @test order2[1]["weights"] == [0, 2, -12]

                                order3 = filter(e -> e["n_sites"] == 3, data)
                                @test length(order3) == 1
                                @test length(order3[1]["bonds"]) == 2
                                @test order3[1]["bonds"][1][end] == 1
                                @test order3[1]["weights"] == [0, 0, 6]
                        end
                end

                @testset "StrongClusterExpansionLattice" begin
                        m_order = 2
                        lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters) + n_unique_sites(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["order"] == 1, data)
                                @test !isempty(order1)
                                @test all(e -> e["n_sites"] == 1, order1)
                                @test length(order1[1]["bonds"]) == 0
                                @test order1[1]["weights"] == [1, -1, 0]

                                order2 = filter(e -> e["order"] == 2, data)
                                @test length(order2) == 1
                                @test length(order2[1]["bonds"]) == 6
                                @test order2[1]["bonds"][1][end] == 1
                                @test order2[1]["weights"] == [0, 0.25, -3.0]

                                order3 = filter(e -> e["order"] == 3, data)
                                @test length(order3) == 1
                                @test length(order3[1]["bonds"]) == 13
                                @test order3[1]["bonds"][1][end] == 1
                                @test order3[1]["weights"] == [0, 0, 1.5]
                        end
                end

                @testset "WeakClusterExpansionLattice" begin
                        m_order = 2
                        lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters) + n_unique_sites(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["order"] == 1, data)
                                @test !isempty(order1)
                                @test all(e -> e["n_sites"] == 1, order1)
                                @test length(order1[1]["bonds"]) == 0
                                @test order1[1]["weights"] == [1, -2, 1]

                                order2 = filter(e -> e["order"] == 2, data)
                                @test length(order2) == 1
                                @test length(order2[1]["bonds"]) == 4
                                @test order2[1]["bonds"][1][end] == 1
                                @test order2[1]["weights"] == [0, 0.5, -2]

                                order3 = filter(e -> e["order"] == 3, data)
                                @test length(order3) == 1
                                @test length(order3[1]["bonds"]) == 8
                                @test order3[1]["bonds"][1][end] == 1
                                @test order3[1]["weights"] == [0, 0, 1]
                        end
                end
        end

        @testset "write_to_json isomorphic + symmetric relation" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, square_uc)

                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)

                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)
                sym_clusters = SymmetricClusterSet(lattice, :Square)
                clusters_from_clusters!(sym_clusters, trans_clusters)

                iso_expansion = Expansion(iso_clusters, lattice)
                summation!(iso_expansion, m_order)
                sym_expansion = Expansion(sym_clusters, lattice)
                summation!(sym_expansion, m_order)

                adj = bond_matrix(lattice)
                colors = Lincege.get_site_colors(lattice)
                canon = Dict{UInt,Tuple}()
                for c in Iterators.flatten((collect(iso_clusters), collect(sym_clusters)))
                        h, perm = Lincege.ghash_with_permutation(iso_clusters.hasher, c.vertices)
                        order = sortperm(collect(perm))
                        v = collect(c.vertices)
                        form = (adj[v, v][order, order], collect(colors[v])[order])
                        if haskey(canon, h)
                                @test canon[h] == form
                        else
                                canon[h] = form
                        end
                end

                mktempdir() do dir
                        path = joinpath(dir, "related.json")
                        write_to_json(iso_expansion, iso_clusters, sym_expansion, sym_clusters, lattice, path)
                        data = JSON.parse(read(path, String))

                        @test length(data) == length(iso_clusters)

                        for entry in data
                                @test haskey(entry, "canonical_bonds")
                                @test haskey(entry, "canonical_site_colors")
                                @test haskey(entry, "related_symmetric")
                        end

                        total_related = sum(length(entry["related_symmetric"]) for entry in data)
                        @test total_related == length(sym_clusters)

                        for entry in data
                                n = entry["n_sites"]
                                for rel in entry["related_symmetric"]
                                        @test length(rel["permuted_coordinates"]) == n
                                        @test haskey(rel, "weights")
                                end
                        end
                end
        end

        @testset "import_from_json round-trip (SiteExpansionLattice)" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)
                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                mktempdir() do dir
                        path = joinpath(dir, "expansion.json")
                        write_to_json(expansion, lattice, path)
                        raw = JSON.parse(read(path, String))
                        imported = import_from_json(path)

                        @test length(imported) == length(iso_clusters)
                        @test length(imported) == length(raw)

                        expansion_hashes = Set(Iterators.flatten(expansion.order_ids))
                        @test Set(c.cluster_hash for c in imported) == expansion_hashes

                        raw_by_hash = Dict(e["cluster_hash"] => e for e in raw)
                        for c in imported
                                key = string(c.cluster_hash)
                                @test haskey(raw_by_hash, key)
                                e = raw_by_hash[key]
                                @test c.order == e["order"]
                                @test c.n_sites == e["n_sites"]
                                @test c.site_colors == Int.(e["site_colors"])
                                @test length(c.bonds) == length(e["bonds"])
                                for (b, eb) in zip(c.bonds, e["bonds"])
                                        @test [b[1], b[2], b[3]] == Int.(eb)
                                end
                                @test isapprox(c.weights, Float64.(e["weights"]))
                        end
                end
        end

        @testset "import_from_json_by_order (SiteExpansionLattice)" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)
                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                mktempdir() do dir
                        path = joinpath(dir, "expansion.json")
                        write_to_json(expansion, lattice, path)

                        by_order = import_from_json_by_order(path)
                        flat = import_from_json(path)

                        for order in eachindex(by_order)
                                @test all(c -> c.order == order, by_order[order])
                        end
                        @test sum(length, by_order) == length(flat)
                end
        end

        @testset "Resummations" begin
                # Alternating harmonic series: partial sums converge to log(2).
                n_orders = 30
                terms = [(-1.0)^(n + 1) / n for n in 1:n_orders]
                energy = [[t] for t in terms]
                specific_heat = [[0.0] for _ in terms]
                entropy = [[0.0] for _ in terms]
                temperatures = [1.0]
                result = NLCEResult(temperatures, energy, specific_heat, entropy, n_orders)

                @testset "euler_resummation" begin
                        resum = euler_resummation(result, 1)
			@test length(resum.energy) == n_orders + 1
                        @test isapprox(abs(resum.energy[end][1]), log(2), atol = 1e-6)
                end

                @testset "wynn_resummation" begin
                        resum = wynn_resummation(result, 6)
			@test length(resum.energy) == n_orders + 1
                        @test isapprox(resum.energy[end][1], log(2), atol = 1e-6)
                end

                @testset "apply_resummations" begin
                        chained = apply_resummations(result, [(:Euler, 1), (:Wynn, 5)])
			@test length(chained.energy) == n_orders + 2
                        @test isapprox(chained.energy[end][1], log(2), atol = 1e-6)
                        @test_throws ErrorException apply_resummations(result, [(:Unknown, 1)])
                end
        end

end # Expansions
