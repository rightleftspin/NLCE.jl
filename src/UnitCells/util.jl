
function label_tiling_units(exp_basis::AbstractVector{<:AbstractMatrix{<:Real}}, primitive_vectors::AbstractMatrix{<:Real})
        pw_dir = pairwise_direction(hcat(exp_basis...))

        labels = find_labels(pw_dir, primitive_vectors)

        partition_labels(labels, exp_basis)
end

function find_labels(pw_dir::AbstractArray{<:Real,3}, primitive_vectors::AbstractMatrix{<:Real})
        n_coords, _, _ = size(pw_dir)
        labels = Vector(1:n_coords)
        prim_vector_math = transpose(primitive_vectors)

        for i in 1:n_coords
                for j in i:n_coords
                        direction = pw_dir[i, j, :]
                        int_coeffs = prim_vector_math \ direction
                        if all(x -> isapprox(x, round(x)), int_coeffs)
                                m = min(labels[i], labels[j])
                                labels[i] = m
                                labels[j] = m
                        end

                end
        end

        labels
end

function partition_labels(labels::AbstractVector{<:Int}, exp_basis::AbstractVector{<:AbstractMatrix{<:Real}})
        partition_sizes = size.(exp_basis, 2)
        partitioned_labels = Vector{Vector{Int}}()
        start_idx = 1
        for partition_size in partition_sizes
                push!(partitioned_labels, labels[start_idx:partition_size+start_idx-1])
                start_idx += partition_size
        end

        partitioned_labels
end


function is_strong_tiling(unit_cell::ExpansionUnitCell)
        all_labels = collect(Iterators.flatten(unit_cell.translation_labels))
        length(all_labels) == length(unique(all_labels))
end