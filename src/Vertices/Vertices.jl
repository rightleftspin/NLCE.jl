"""
    AbstractVertices{V}

Abstract base type for a set of vertices of element type `V`.

Subtypes must implement:
- `vertices(vertex_set)` - return the underlying iterable of vertex indices
- `Base.collect(vertex_set)` - return a sorted `Vector{V}`
- `Base.sort(vertex_set)` - return a sorted copy
- `Base.intersect(vertex_set1, vertex_set2)` - set intersection
- `Base.setdiff(vertex_set1, vertex_set2)` - set difference
- `Base.union(vertex_set1, vertex_set2)` - set union (two-argument form)
- `Base.in(v, vertex_set)` - membership test
- `Base.eltype(vertex_set)` - element type `V`
"""
abstract type AbstractVertices{V} end

vertices(vertex_set::AbstractVertices) = _NI("vertices")
Base.collect(vertex_set::AbstractVertices) = _NI("Base.collect")
Base.sort(vertex_set::AbstractVertices) = _NI("Base.sort")
Base.intersect(vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = _NI("Base.intersect")
Base.setdiff(vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = _NI("Base.setdiff")
Base.union(vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = _NI("Base.union")
Base.in(v, vertex_set::AbstractVertices) = _NI("Base.in")
Base.eltype(vertex_set::AbstractVertices) = _NI("Base.eltype")

Base.length(vertex_set::AbstractVertices) = length(vertices(vertex_set))
Base.isempty(vertex_set::AbstractVertices) = isempty(vertices(vertex_set))

Base.isequal(vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = vertex_set1 == vertex_set2
Base.:(==)(vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = (collect(vertex_set1) == collect(vertex_set2))
Base.hash(vertex_set::AbstractVertices, h::UInt) = hash(collect(sort(vertex_set)), h)
Base.contains(vertex_set::AbstractVertices, v) = v in vertex_set
Base.haskey(vertex_set::AbstractVertices, v) = v in vertex_set

function Base.union(vertex_set::AbstractVertices, itr)
        vertex_set_temp = vertex_set
        for x in itr
                vertex_set_temp = union(vertex_set_temp, x)
        end
        vertex_set_temp
end

Base.getindex(vec::Vector, vertex_set::AbstractVertices) = vec[collect(vertex_set)]
Base.getindex(mat::Matrix, vertex_set1::AbstractVertices, vertex_set2::AbstractVertices) = mat[collect(vertex_set1), collect(vertex_set2)]
Base.getindex(mat::Matrix, ::Colon, vertex_set::AbstractVertices) = mat[:, collect(vertex_set)]
Base.getindex(mat::Matrix, vertex_set::AbstractVertices, ::Colon) = mat[collect(vertex_set), :]

Base.iterate(vertex_set::AbstractVertices) = iterate(vertices(vertex_set))
Base.iterate(vertex_set::AbstractVertices, state) = iterate(vertices(vertex_set), state)

Base.show(io::IO, vertex_set::AbstractVertices) = print(io, "Vertices: ", collect(vertex_set))

include("TaggedVertices.jl")
