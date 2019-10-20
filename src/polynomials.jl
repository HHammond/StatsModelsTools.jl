"""Polynomial regression interface

Taken from https://juliastats.github.io/StatsModels.jl/stable/internals/#extending-1
"""

module Polynomials

using StatsModels
using StatsBase

poly(x, n) = x.^n

const POLY_CONTEXT = Any

struct PolyTerm{T,D} <: AbstractTerm
    term::T
    deg::D
end

Base.show(io::IO, p::PolyTerm) = print(io, "poly($(p.term), $(p.deg))")

poly(t::Symbol, d::Int) = PolyTerm(term(t), term(d))

function StatsModels.apply_schema(t::FunctionTerm{typeof(poly)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:POLY_CONTEXT})
    apply_schema(PolyTerm(t.args_parsed...), sch, Mod)
end

function StatsModels.apply_schema(t::PolyTerm,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:POLY_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("PolyTerm only works with continuous terms (got $term)"))
    isa(t.deg, ConstantTerm) ||
        throw(ArgumentError("PolyTerm degree must be a number (got $t.deg)"))
    PolyTerm(term, t.deg.n)
end

function StatsModels.modelcols(p::PolyTerm, d::NamedTuple)
    col = modelcols(p.term, d)
    reduce(hcat, [col.^n for n in 1:p.deg])
end

StatsModels.terms(p::PolyTerm) = terms(p.term)
StatsModels.termvars(p::PolyTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::PolyTerm) = p.deg

StatsBase.coefnames(p::PolyTerm) = coefnames(p.term) .* "^" .* string.(1:p.deg)

export poly
end
