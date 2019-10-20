"""Basic interface for spline functions
"""

using StatsModels
using StatsBase

const SplineTermContext = Any

abstract type AbstractSplineTerm <: AbstractTerm end

shorthand(::Type{<:AbstractSplineTerm}) = "spline"
spline_term(v::AbstractSplineTerm) = throw("term not implemented")
spline_knots(v::AbstractSplineTerm) = throw("knots not implemented")
spline_deg(v::AbstractSplineTerm) = throw("deg not implemented")

function Base.show(io::IO, t::T) where T<:AbstractSplineTerm
    print(io, "$(shorthand(T))($(spline_term(t)), $(spline_knots(t)))")
end

function StatsBase.coefnames(t::T) where T <: AbstractSplineTerm
    c = coefnames(spline_term(t))
    k = 1:width(t)

    map(k) do knot
        "$(shorthand(T))($(c), $(t.knots.n))[$knot]"
    end
end
