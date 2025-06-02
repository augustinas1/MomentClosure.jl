abstract type MomentEquations end

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
    """Dictionary mapping species of the associated `ReactionSystem` to their index in the moment equations."""
    smap::Dict
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
    """Dictionary mapping species of the associated `ReactionSystem` to their index in the moment equations."""
    smap::Dict
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

"""
```julia
SciMLBase.ODEProblem(sys::MomentEquations, u0, tspan, p=NullParameters(); 
                     use_deterministic_IC::Bool=true, kwargs...)
```

Generates a [`SciMLBase.ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem) from `MomentEquations`.

Note that the initial condition `u0` can be defined as follows:
- If `u0` contains the initial molecule numbers (either as a vector or a mapping) and `use_deterministic_IC` is set to `true`, 
  `u0` will be mapped to the corresponding raw/central moments under the assumption of deterministic 
  initial conditions (using the function [`deterministic_IC`](@ref)). This is the default use case.
- To handle a more complicated initial condition, `u0` can be defined with the initial values of each moment directly, 
  but `use_deterministic_IC` must be set to `false` to ensure correctness.
"""
function SciMLBase.ODEProblem(sys::MomentEquations, u0, tspan, p=NullParameters(); 
                              use_deterministic_IC::Bool=true, kwargs...)
    u0map = use_deterministic_IC ? deterministic_IC(u0, sys) : u0
    ODEProblem(complete(get_odes(sys)), u0map, tspan, p; kwargs...)
end

function Base.nameof(eqs::MomentEquations)
    nameof(get_odes(eqs))
end

# Basic `MomentEquations`-specific accessors
get_odes(sys::MomentEquations) = getfield(sys, :odes)
get_closure(sys::ClosedMomentEquations) = getfield(sys, :closure)
Catalyst.speciesmap(sys::MomentEquations) = getfield(sys, :smap)
Catalyst.speciesmap(sys::ClosedMomentEquations) = speciesmap(getfield(sys, :open_eqs))
ModelingToolkit.get_iv(sys::MomentEquations) = get_iv(get_odes(sys))
ModelingToolkit.get_eqs(sys::MomentEquations) = get_eqs(get_odes(sys))
ModelingToolkit.unknowns(sys::MomentEquations) = unknowns(get_odes(sys))
ModelingToolkit.get_ps(sys::MomentEquations) = get_ps(get_odes(sys))