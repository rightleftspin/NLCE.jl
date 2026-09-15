# Utility functions that will be used in multiple modules
function pairwise_direction(coordinates::AbstractMatrix{<:Real})
        T = eltype(coordinates)
        dim, n_coords = size(coordinates)
        pw_dir = Array{T}(undef, n_coords, n_coords, dim)
        @inbounds for i in 1:n_coords
                for j in 1:n_coords
                        for d in 1:dim
                                pw_dir[i, j, d] = coordinates[d, j] - coordinates[d, i]
                        end
                end
        end
        pw_dir
end

const all_lattice_symmetries = Dict(
        :Square => Vector{Matrix{Float64}}([
                [1 0; 0 1],
                [0 1; -1 0],
                [-1 0; 0 -1],
                [0 -1; 1 0],
                [0 1; 1 0],
                [0 -1; -1 0],
                [1 0; 0 -1],
                [-1 0; 0 1],
        ]),
        :Shastry_Sutherland => Vector{Matrix{Float64}}([
                [1 0; 0 1],
                [1 0; 0 -1],
                [-1 0; 0 1],
        ]),
)

# ---------------------------------------------------------------------------
# Example unit cells
# ---------------------------------------------------------------------------

"""
    square_unit_cell

Square-lattice unit cell (site expansion).
"""
const square_unit_cell = UnitCell(
        [[0.0, 0.0]],
        [[1.0, 0.0], [0.0, 1.0]],
        [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)],
        [1],
)

"""
    kagome_unit_cell

Kagome-lattice unit cell (site expansion).
"""
const kagome_unit_cell = UnitCell(
        [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]],
        [[2.0, 0.0], [1.0, sqrt(3)]],
        [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
                Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)],
        [1, 1, 1],
)

"""
    pyrochlore_expansion_unit_cell

Pyrochlore unit-cell expansion (`ExpansionUnitCell`, strong cluster expansion).
"""
const pyrochlore_expansion_unit_cell = ExpansionUnitCell(
        [
                [
                        [1 / 2, 1 / 2, 1 / 2],
                        [1 / 2, -1 / 2, -1 / 2],
                        [-1 / 2, 1 / 2, -1 / 2],
                        [-1 / 2, -1 / 2, 1 / 2],
                ],
        ],
        [[2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]],
        [
                # Inside a unit cell
                ExpansionBond([1, 1], [1, 2], [0, 0, 0], 1),
                ExpansionBond([1, 1], [1, 3], [0, 0, 0], 1),
                ExpansionBond([1, 1], [1, 4], [0, 0, 0], 1),
                ExpansionBond([1, 2], [1, 3], [0, 0, 0], 1),
                ExpansionBond([1, 2], [1, 4], [0, 0, 0], 1),
                ExpansionBond([1, 3], [1, 4], [0, 0, 0], 1),
                # Between unit cells
                ExpansionBond([1, 1], [1, 2], [0, 0, 1], 1),
                ExpansionBond([1, 1], [1, 3], [0, 1, 0], 1),
                ExpansionBond([1, 1], [1, 4], [1, 0, 0], 1),
                ExpansionBond([1, 2], [1, 3], [0, 1, -1], 1),
                ExpansionBond([1, 2], [1, 4], [1, 0, -1], 1),
                ExpansionBond([1, 3], [1, 4], [1, -1, 0], 1),
        ],
        [
                Bond(1, 1, [1, 0, 0], 1),
                Bond(1, 1, [0, 1, 0], 1),
                Bond(1, 1, [0, 0, 1], 1),
                Bond(1, 1, [0, 1, -1], 1),
                Bond(1, 1, [1, 0, -1], 1),
                Bond(1, 1, [1, -1, 0], 1),
        ],
        [[1, 1, 1, 1]],
)

"""
    square_cluster_expansion_unit_cell

Square-cluster unit-cell expansion (`ExpansionUnitCell`, weak cluster expansion).
"""
const square_cluster_expansion_unit_cell = ExpansionUnitCell(
        [[[-1 / 2, -1 / 2], [-1 / 2, 1 / 2], [1 / 2, -1 / 2], [1 / 2, 1 / 2]]],
        [[1.0, 1.0], [1.0, -1.0]],
        [
                ExpansionBond([1, 1], [1, 2], [0, 0], 1),
                ExpansionBond([1, 1], [1, 3], [0, 0], 1),
                ExpansionBond([1, 2], [1, 4], [0, 0], 1),
                ExpansionBond([1, 3], [1, 4], [0, 0], 1),
        ],
        [
                Bond(1, 1, [1, 0], 1),
                Bond(1, 1, [0, 1], 1),
        ],
        [[1, 1, 1, 1]],
)
