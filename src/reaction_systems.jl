"""
    propensities(rn::ReactionSystem; combinatoric_ratelaws=true)
Return a vector of propensity functions of all reactions in the given [`ReactionSystem`]
(https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.ReactionSystem).

Notes:
- `combinatoric_ratelaws=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`Catalyst.jumpratelaw`]
  (https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.jumpratelaw).
"""
function propensities(rn::ReactionSystem; combinatoric_ratelaws::Bool=true)
    simplify.(jumpratelaw.(reactions(rn); combinatoric_ratelaw=combinatoric_ratelaws))
end

"""
    get_stoichiometry(rn::ReactionSystem, smap::AbstractDict)
Return the net stoichiometry matrix using the specified mapping of species to their indices.

Notes:
- This is a modification of [`Catalyst.netstoichmat`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.netstoichmat) 
  that is used internally to deal with reactions involving symbolic stoichiometry coefficients.
- The function also allows custom `smap`, so it is not limited to the default 
  [`Catalyst.speciesmap`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.speciesmap) ordering.
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