function generate_coordinates(dim_specs::NTuple{N,Union{Int,UnitRange{Int}}}, num_basis_elements::Int) where {N}
        ranges = [s isa Int ? (-s:s) : s for s in dim_specs]
        primitive = hcat([collect(t) for t in Iterators.product(ranges...)]...)
        hcat([vcat(primitive, fill(i, 1, size(primitive, 2))) for i in 1:num_basis_elements]...)
end

generate_coordinates(max_order::Int, num_basis_elements::Int, dimension::Int) =
        generate_coordinates(ntuple(_ -> max_order, dimension), num_basis_elements)

function build_lattice_coordinates(max_order::Int, expansion_unit_cell::ExpansionUnitCell;
        coord_gen=(specs, n) -> generate_coordinates(specs, n))
        d = dimension(expansion_unit_cell)
        blocks = Matrix{Int}[]
        for i in 1:length(basis_size(expansion_unit_cell))
                coords = coord_gen(ntuple(_ -> max_order, d), basis_size(expansion_unit_cell)[i])
                push!(blocks, vcat(coords[1:d, :], ones(Int, size(coords, 2))' * i, coords[end, :]'))
        end
        hcat(blocks...)
end


function generate_finite_coordinates(dims::NTuple{N,Int}, num_basis::Int) where {N}
        generate_coordinates(map(d -> 0:d-1, dims), num_basis)
end

function build_lattice_coordinates_finite(dims::NTuple{N,Int}, expansion_unit_cell::ExpansionUnitCell) where {N}
        build_lattice_coordinates(0, expansion_unit_cell;
                coord_gen=(_, n) -> generate_finite_coordinates(dims, n))
end

function generate_coord_index(coordinates::Matrix{Int})
        coord_index = Dict{Vector{Int},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[col] = i
        end
        coord_index
end

function generate_coord_index(coordinates::Matrix{Float64})
        coord_index = Dict{Vector{Float64},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[round.(col, digits=6)] = i
        end
        coord_index
end

function _fill_adj_matrix!(adj_matrix, coordinates, coord_index, bonds, site_matcher, neighbor_fn=neighbor_site)
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in bonds
                        if site_matcher(bond, col)
                                ni = get(coord_index, neighbor_fn(bond, col), nothing)
                                if ni !== nothing
                                        adj_matrix[index, ni] = bond.bond_type
                                        adj_matrix[ni, index] = bond.bond_type
                                end
                        end
                end
        end
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell; neighbor_fn=neighbor_site)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        _fill_adj_matrix!(adj_matrix, coordinates, coord_index, unit_cell.bonds, (bond, col) -> bond.site1 == col[end], neighbor_fn)
        adj_matrix
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell; neighbor_fn=neighbor_site)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        _fill_adj_matrix!(adj_matrix, coordinates, coord_index, unit_cell.bonds, (bond, col) -> bond.site1 == col[end-1:end], neighbor_fn)
        adj_matrix
end

function generate_adj_matrix_weak(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell, unique_inds::Vector{Int}; neighbor_fn=neighbor_site)
        real_coords = shift_unit_cell(unit_cell, coordinates)[:, unique_inds]
        coord_index = generate_coord_index(real_coords)
        adj_matrix = zeros(Int, length(unique_inds), length(unique_inds))
        for coord in eachcol(coordinates)
                trans_ind = get(coord_index, round.(shift_unit_cell(unit_cell, coord), digits=6), nothing)
                isnothing(trans_ind) && continue
                for bond in unit_cell.bonds
                        if bond.site1 == coord[end-1:end]
                                neighbor_coord = neighbor_fn(bond, coord)
                                ni = get(coord_index, round.(shift_unit_cell(unit_cell, neighbor_coord), digits=6), nothing)
                                if ni !== nothing
                                        adj_matrix[trans_ind, ni] = bond.bond_type
                                        adj_matrix[ni, trans_ind] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function _generate_neighbor_list(coordinates::AbstractMatrix{Int}, bonds, ::Type{V}, neighbor_fn=neighbor_site) where {V<:AbstractVertices}
        coord_index = generate_coord_index(coordinates)
        neighbor_list = fill(V(), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in bonds
                        if bond.site1 == col[end]
                                neighbor_index = get(coord_index, neighbor_fn(bond, col), nothing)
                                if neighbor_index !== nothing
                                        neighbor_list[index] = union(neighbor_list[index], V(neighbor_index))
                                        neighbor_list[neighbor_index] = union(neighbor_list[neighbor_index], V(index))
                                end
                        end
                end
        end

        neighbor_list
end

generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell; neighbor_fn=neighbor_site) =
        _generate_neighbor_list(coordinates, unit_cell.bonds, LatticeVertices{Int}, neighbor_fn)
generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell; neighbor_fn=neighbor_site) =
        _generate_neighbor_list(coordinates, unit_cell.expansion_bonds, ExpansionVertices{Int}, neighbor_fn)

find_centers(coordinates::AbstractMatrix{Int}) = findall(col -> all(==(0), col[1:end-1]), eachcol(coordinates))

function generate_strong_connections(expansion_coordinates::AbstractMatrix{Int}, lattice_coordinates::AbstractMatrix{Int})

        connections_vec = fill(LatticeVertices{Int}(), size(expansion_coordinates, 2))
        lattice_slice = @view lattice_coordinates[1:end-1, :]
        for (i, coord) in enumerate(eachcol(expansion_coordinates))
                connection = findall(==(coord), eachcol(lattice_slice))
                connections_vec[i] = LatticeVertices(collect(connection))
        end
        connections_vec
end

# Source - https://stackoverflow.com/a/50900113
# Posted by Bogumil Kaminski
# Retrieved 2026-04-02, License - CC BY-SA 4.0
# Modified by adding a round function to deal
# with floating point numbers
function uniqueidx(x::AbstractArray{T}) where {T}
        uniqueset = Set{T}()
        ex = eachindex(x)
        idxs = Vector{eltype(ex)}()
        for i in ex
                xi = round.(x[i], digits=6)
                if !(xi in uniqueset)
                        push!(idxs, i)
                        push!(uniqueset, xi)
                end
        end
        idxs
end

function generate_weak_connections(expansion_coordinates::AbstractMatrix{Int}, lattice_coords::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        connections_vec = fill(LatticeVertices{Int}(), size(expansion_coordinates, 2))

        real_coords = shift_unit_cell(unit_cell, lattice_coords)
        unique_inds = uniqueidx(Vector{Vector{Float64}}(eachcol(real_coords)))

        real_coords = real_coords[:, unique_inds]
        reduced_lattice_coordinates = lattice_coords[:, unique_inds]
        reverse_connections_vec = fill(ExpansionVertices{Int}(), size(reduced_lattice_coordinates, 2))

        real_coord_index = Dict(round.(col, digits=6) => i for (i, col) in enumerate(eachcol(real_coords)))
        for (i, coord) in enumerate(eachcol(expansion_coordinates))
                connection = Int[]
                for lcoord in eachcol(lattice_coordinates(unit_cell, coord))
                        con = get(real_coord_index, round.(lcoord, digits=6), nothing)
                        if !isnothing(con)
                                push!(connection, con)
                        end
                end
                connections_vec[i] = LatticeVertices(connection)
                for lv in connection
                        reverse_connections_vec[lv] = union(reverse_connections_vec[lv], ExpansionVertices(i))
                end
        end

        masking_matrix = zeros(Int, size(reduced_lattice_coordinates, 2), size(reduced_lattice_coordinates, 2))
        for (ev_idx, lattice_vertices) in enumerate(connections_vec)
                for lv1 in lattice_vertices, lv2 in lattice_vertices
                        lv1 != lv2 && (masking_matrix[lv1, lv2] = ev_idx)
                end
        end

        return connections_vec, reverse_connections_vec, masking_matrix, unique_inds
end




