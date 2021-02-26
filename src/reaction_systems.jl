# System of chemical reactions similar in structure to ModelingToolkit.jl ReactionSystem (parts of code borrowed from there)
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

# Example
```julia
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

# extending a number of ReactionSystem functions to support ReactionSystemMod

#import Catalyst: species, speciesmap, numspecies, numreactions, numparams
function species(rn::ReactionSystemMod)
    rn.states
end

function params(rn::ReactionSystemMod)
    rn.ps
end

function speciesmap(rn::ReactionSystemMod)
    Dict(S => i for (i,S) in enumerate(species(rn)))
end

function paramsmap(rn::ReactionSystemMod)
    Dict(p => i for (i,p) in enumerate(params(rn)))
end

function numspecies(rn::ReactionSystemMod)
    length(rn.states)
end

function numreactions(rn::ReactionSystemMod)
    length(rn.a)
end

function numparams(rn::ReactionSystemMod)
    length(rn.ps)
end

function get_S_mat(rn::Union{ReactionSystem, ReactionSystemMod})
    if typeof(rn) == ReactionSystem
        prodstoichmat(rn)' - substoichmat(rn)'
    else
        rn.S
    end
end

function propensities(rn::Union{ReactionSystem, ReactionSystemMod}; combinatoric_ratelaw=true)
    if typeof(rn) == ReactionSystem
        [simplify(jumpratelaw(rx, combinatoric_ratelaw=combinatoric_ratelaw)) for rx in reactions(rn)]
    else
        rn.a
    end
end
