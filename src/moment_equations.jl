abstract type MomentEquations end

# TODO: implement getters

"""
$(TYPEDEF)

Raw moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct RawMomentEquations <: MomentEquations
    """[`ModelingToolkit.ODESystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/#ModelingToolkit.ODESystem)
    consisting of the time-evolution equations of raw moments."""
    odes::ODESystem
    """Symbolic variables defining the raw moments."""
    μ::OrderedDict
    """Number of species within the system."""
    N::Int
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector
    """Vector of all index combinations up to `m_order`."""
    iter_m::Vector
    """Vector of all index combinations of order greater than `m_order`
    up to `q_order`."""
    iter_q::Vector
    """Vector of index combinations of order 1."""
    iter_1::Vector
end

"""
$(TYPEDEF)

Central moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct CentralMomentEquations <: MomentEquations
    """[`ModelingToolkit.ODESystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/#ModelingToolkit.ODESystem)
    consisting of the time-evolution equations of central moments."""
    odes::ODESystem
    """Symbolic variables defining the means."""
    μ::OrderedDict
    """Symbolic variables defining the central moments."""
    M::Dict
    """Number of species within the system."""
    N::Int
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector
    """Vector of all index combinations up to `m_order`."""
    iter_m::Vector
    """Vector of all index combinations of order greater than `m_order`
    up to `q_order`."""
    iter_q::Vector
    """Vector of index combinations of order 1."""
    iter_1::Vector
end

"""
$(TYPEDEF)

Closed moment equations and the corresponding closure functions.

# Fields
$(FIELDS)
"""
struct ClosedMomentEquations <: MomentEquations
    """[`ModelingToolkit.ODESystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/#ModelingToolkit.ODESystem)
    consisting of the time-evolution equations of *closed* moments."""
    odes::ODESystem
    """Dictionary of moment closure functions for each higher order moment."""
    closure::OrderedDict
    """Original raw or central moment equations (before closure was applied)."""
    open_eqs::MomentEquations
end

# a basic wrapper
function SciMLBase.ODEProblem(eqs::MomentEquations, args...; kwargs...)
    ODEProblem(complete(eqs.odes), args...; kwargs...)
end

function Base.nameof(eqs::MomentEquations)
    nameof(eqs.odes)
end
