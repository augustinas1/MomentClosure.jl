"""
    propensities(rn::Union{ReactionSystem, ReactionSystemMod}; combinatoric_ratelaw=true)
Return a vector of propensity functions of all reactions in the given [`ReactionSystem`]
(https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem).

Notes:
- `combinatoric_ratelaw=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
"""
function propensities(rn::ReactionSystem; combinatoric_ratelaw=true)
    simplify.(jumpratelaw.(reactions(rn); combinatoric_ratelaw))
end

"""
    get_stoichiometry(rn::ReactionSystem, smap::AbstractDict)
Return the net stoichiometry matrix using the specified mapping of species to their indices.

Notes:
- This is a modification of [`Catalyst.netstoichmat`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.netstoichmat) 
  that is used internally to deal with reactions involving symbolic stoichiometry coefficients.
- The function also allows custom `smap`, so it is not limited to the default 
  [`Catalyst.speciesmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.speciesmap) ordering.
- TODO: remove once [this Catalyst issue](https://github.com/SciML/Catalyst.jl/issues/489)
  is resolved.
"""
function get_stoichiometry(rn::ReactionSystem, smap::AbstractDict)
    nmat = Matrix(undef, numspecies(rn), numreactions(rn))
    fill!(nmat, 0)
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
            nmat[smap[spec],k] = coef
        end
    end
    nmat
end