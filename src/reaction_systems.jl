# System of chemical reactions similar in structure to ModelingToolkit.jl ReactionSystem
# allow stochastic variables (sampled from the geometric distribution) as stoichometry matrix coefficients (hence Mod)
# Note that this construction is very limited compared to ReactionSystem and is not integrated within the broader framework of ModelingToolkit.jl
# parts of code adapted from
# https://github.com/SciML/ModelingToolkit.jl/blob/master/src/systems/reaction/reactionsystem.jl
"""
$(TYPEDEF)

A system of chemical reactions. The structure is heavily based on
`ModelingToolkit.ReactionSystem` and the relevant fields are left
identical for the sake of familiar users.

# Fields
$(FIELDS)

# Notes
- The molecule numbers of each chemical species must be explicitly defined as
  [`Symbolics.@variables`](@ref) and reaction constants as [`ModelingToolkit.@parameters`](@ref).
  The time, ``t``, must also be initialised as a `ModelingToolkit.parameter`
  and be used to indicate the time-dependence of species' variables.
- To specify that a certain reaction product is geometrically distributed, the
  corresponding stoichiometric matrix element must be a `ModelingToolkit.parameter`
  denoting the mean of the distribution.

# Example
```julia
using MomentClosure

@parameters t, c₁, c₂, c₃, c₄, Ω
@variables X(t), Y(t)

# stoichiometric matrix
S_mat = [ 1 -1  1 -1;
         -1  1  0  0]

# propensity functions
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]

rn = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)
```
"""
struct ReactionSystemMod
    """Independent variable (usually time)"""
    iv::Any
    """Dependent (state) variables representing amount of each species."""
    states::Vector
    """Parameter variables."""
    ps::Vector
    """Propensity functions of all reactions"""
    a::Vector
    """Stoichiometric matrix"""
    S::Matrix

    function ReactionSystemMod(iv, states, ps, a, S)
        if size(S)[1] != size(states)[1]
            error("Inconsistent stoichiometric matrix dimensions and number of species")
        elseif size(S)[2] != size(a)[1]
            error("Inconsistent stoichiometric matrix dimensions and number of reactions")
        else
            if (typeof(iv) != Num)
                error("the indepent variable is not Num type (possibly breaking update of Catalyst/ModelingToolkit)")
            else
            # initially using @parameters and @variables macros we have defined Num type symbols
            # which have to be converted to SymbolicUtils supported types for further manipulations
                new(value(iv), value.(states), value.(ps), simplify.(value.(a)), value.(S))
            end
        end
    end

end

# Implementing a number of basic functions allowing easy access to reaction network properties.
# The functions are identical to the Catalyst ones (full credit to the developers) in
# order to keep the APIs consistent
# TODO: use getters

"""
    species(rn::ReactionSystemMod)
Given a [`ReactionSystemMod`](@ref), return a vector of species variables expressed as
[`Term{Real}`](https://mtk.sciml.ai/stable/IR/#Types-1). Extension of Catalyst's
[`species`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.species) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function species(rn::ReactionSystemMod)
    rn.states
end

"""
    params(rn::ReactionSystemMod)
Given a [`ReactionSystemMod`](@ref), return a vector of parameter variables
expressed as `Sym{ModelingToolkit.Parameter{Real}`. Extension of  Catalyst's
[`params`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.params) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function params(rn::ReactionSystemMod)
    rn.ps
end

"""
    speciesmap(rn::ReactionSystemMod)
Given a [`ReactionSystemMod`](@ref), return a Dictionary mapping from species to
species indices. Extension of Catalyst's
[`speciesmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.speciesmap) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function speciesmap(rn::ReactionSystemMod)
    Dict(S => i for (i,S) in enumerate(species(rn)))
end

"""
    paramsmap(rn::ReactionSystemMod)
Given a [`ReactionSystemMod`](@ref), return a Dictionary mapping from parameters to
parameter indices.  Extension of Catalyst's
[`paramsmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.paramsmap) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function paramsmap(rn::ReactionSystemMod)
    Dict(p => i for (i,p) in enumerate(params(rn)))
end

"""
    numspecies(rn::ReactionSystemMod)
Return the number of species within the given [`ReactionSystemMod`](@ref). Extension of Catalyst's
[`numspecies`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.numspecies) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function numspecies(rn::ReactionSystemMod)
    length(rn.states)
end

"""
    numreactions(rn::ReactionSystemMod)
Return the number of reactions within the given [`ReactionSystemMod`](@ref). Extension of Catalyst's
[`numreactions`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.numreactions) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function numreactions(rn::ReactionSystemMod)
    length(rn.a)
end

"""
    numparams(rn::ReactionSystemMod)
Return the number of parameters within the given [`ReactionSystemMod`](@ref). Extension of Catalyst's
[`numparams`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.numparams) used with a
[`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).
"""
function numparams(rn::ReactionSystemMod)
    length(rn.ps)
end

"""
    get_S_mat(rn::Union{ReactionSystem, ReactionSystemMod})
Return the (net) stoichiometric matrix of the given [`ReactionSystem`]
(https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
or [`ReactionSystemMod`](@ref).
"""
function get_S_mat(rn::Union{ReactionSystem, ReactionSystemMod})
    if typeof(rn) == ReactionSystem
        prodstoichmat(rn)' - substoichmat(rn)'
    else
        rn.S
    end
end

"""
    propensities(rn::Union{ReactionSystem, ReactionSystemMod}; combinatoric_ratelaw=true)
Return a vector of propensity functions of all reactions in the given [`ReactionSystem`]
(https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
or [`ReactionSystemMod`](@ref).

Notes:
- `combinatoric_ratelaw=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
  *Note* that this field is irrelevant using `ReactionSystemMod` as then the
  propensities are defined directly by the user.
"""
function propensities(rn::Union{ReactionSystem, ReactionSystemMod}; combinatoric_ratelaw=true)
    if typeof(rn) == ReactionSystem
        simplify.(jumpratelaw.(reactions(rn), combinatoric_ratelaw=combinatoric_ratelaw))
    else
        rn.a
    end
end
