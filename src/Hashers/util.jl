function unique_direction_indices(pw_dir::AbstractArray{<:Real,3}, bond_matrix::AbstractMatrix{Int})
        n_coords = size(pw_dir, 1)
        dir_to_index = Dict{Vector{Float64},Int}()
        index_matrix = zeros(Int, n_coords, n_coords)

        for i in 1:n_coords, j in 1:n_coords
                if i != j && bond_matrix[i, j] != 0
                        dir = round.(pw_dir[i, j, :], digits=6)
                        if !haskey(dir_to_index, dir)
                                dir_to_index[dir] = length(dir_to_index) + 1
                        end
                        index_matrix[i, j] = dir_to_index[dir]
                end
        end

        return index_matrix, length(dir_to_index)
end

function get_permutations(coords::Matrix{Float64}, syms::Vector{Matrix{Float64}})
        n_sites = size(coords, 2)

        # Center coordinates for symmetry operations
        centered = coords .- sum(coords, dims=2) ./ n_sites
        coord_index = Dict(round.(centered[:, i], digits=6) => i for i in 1:n_sites)

        permutations = Vector{Vector{Int64}}()
        for sym in syms
                perm = zeros(Int64, n_sites)
                transformed = sym * centered
                for i in 1:n_sites
                        j = get(coord_index, round.(transformed[:, i], digits=6), nothing)
                        !isnothing(j) && (perm[i] = j)
                end
                push!(permutations, perm)
        end

        permutations
end
