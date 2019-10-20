"""
Regression Modelling Splines from Frank Harrell's RMS book

These methods use truncated power bases to implement regression splines.


Frank E. Harrell, Jr.. 2006. Regression Modeling Strategies. Springer-Verlag, Berlin, Heidelberg.
"""
module RMSSplines

__precompile__()

include("spline_interface.jl")

using StatsModels
using StatsBase


"""Compute Truncated power basis"""
@inline U₊(v) = v > 0 ? v : zero(v)
@inline P(v, k) = v > 0 ? v ^ k : zero(v)
@inline P₃(v) = P(v, 3)

function rms_linear_spline(x::AbstractVector, knots::AbstractVector)
    k = length(knots)

    M = Array{Float64}(undef, k, length(x))

    for i = 1:k
        M[i, :] = U₊.(x .- knots[i])
    end

    M
end

@inline _rcs(x::Missing, j::Int, t::AbstractVector) = missing
@inline _rcs(x::Int, j::Int, t::AbstractVector) = _rcs(float(x), j, t)

"""Compute Truncated Power Basis for Natural/Restricted Cubic Spline """
@inline function _rcs(x, j::Int, t::AbstractVector)
    k = length(t)

    d = (t[k] - t[k-1])
    τ = (t[k] - t[1]) ^ 2

    (
     P₃(x - t[j])
     - P₃(x - t[k - 1]) * (t[k] - t[j]) / d
     + P₃(x - t[k]) * (t[k-1] - t[j]) / d
    ) / τ
end

const RealOrMissing = Union{<:Real, Missing}

@inline function rms_restricted_cubic_spline!(out::AbstractVector, x::Missing, knots::AbstractVector)
    out .= missing
    out
end

@inline function rms_restricted_cubic_spline!(out::AbstractVector, x::Real, knots::AbstractVector)
    k = length(knots)

    out[1] = x
    @simd for j in 2:k-1
        out[j] = _rcs(x, j - 1, knots)
    end

    out
end

function rms_restricted_cubic_spline(x::AbstractVector{<:Real}, knots::AbstractVector{<:Real})
    k = length(knots)
    M = Array{Float64}(undef, length(x), k-1)

    for i in eachindex(x)
        row = view(M, i, :)
        rms_restricted_cubic_spline!(row, x[i], knots)
    end

    M
end

function rms_restricted_cubic_spline(x::AbstractVector{<:RealOrMissing}, knots::AbstractVector{<:Real})
    k = length(knots)

    R = eltype(x) <: Real ? Float64 : Union{Float64, Missing}
    M = Array{R}(undef, length(x), k-1)

    for i in eachindex(x)
        row = view(M, i, :)
        rms_restricted_cubic_spline!(row, x[i], knots)
    end

    M
end

"""RMS quantiles for regression spline knots"""
function rms_quantiles_for_knots(k::Int)
    if k <= 3
        outer = 0.1
    elseif k <= 6
        outer = 0.05
    else
        outer = 0.025
    end

    collect(range(outer, 1-outer, length=k))
end

function rms_generate_knots(x::AbstractVector, knots=3; unique_xs=false)
    if length(x) < knots
        throw("Not enough data to generate knots")
    end

    q = rms_quantiles_for_knots(knots)

    if unique_xs
        quantile(x |> skipmissing |> collect, q)
    else
        quantile(x |> skipmissing |> unique |> collect, q)
    end
end

rcs(x::AbstractVector, knots::Int) = rcs(x, rms_generate_knots(x, knots))
rcs(x::AbstractVector, knots::AbstractVector) = rms_restricted_cubic_spline(x, knots)

################################################################################
"""Define Formula interface for Restricted Cubic Splines"""

struct RCSTerm{T, K} <: AbstractSplineTerm
    term::T
    knots::K
end

shorthand(::Type{<:RCSTerm}) = "rcs"
spline_term(v::RCSTerm) = v.term
spline_knots(v::RCSTerm) = v.knots
spline_deg(v::RCSTerm) = 3

rcs(x::Symbol, knots::Int=3) = RCSTerm(term(x), term(knots))
rcs(t::AbstractTerm, knots::Int=3) = RCSTerm(t, term(knots))

function StatsModels.apply_schema(t::FunctionTerm{typeof(rcs)}, sch::StatsModels.Schema, Mod::Type{<:SplineTermContext})
    apply_schema(RCSTerm(t.args_parsed...), sch, Mod)
end

function StatsModels.apply_schema(t::RCSTerm, sch::StatsModels.Schema, Mod::Type{<:SplineTermContext})
    term = apply_schema(t.term, sch, Mod)

    isa(t.knots, ConstantTerm) || throw("Knots must be a constant term")
    t.knots.n > 2 || throw("Restricted cubic spline requires at least 3 knots")

    RCSTerm(term, t.knots)
end

function StatsModels.modelcols(t::RCSTerm, d::NamedTuple)
    col = modelcols(t.term, d)
    rcs(col, t.knots.n)
end

StatsModels.width(t::RCSTerm) = t.knots.n - 1

export rcs, rms_generate_knots, rms_quantiles_for_knots

end
