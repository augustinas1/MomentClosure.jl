"""
      linear_mapping_approximation(rn_nonlinear::T, rn_linear::T, binary_vars::Array{Int,1}=Int[], m_order::Int=0;
                                   combinatoric_ratelaw = true) where T <: ReactionSystem

Given a *nonlinear* [`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
and an equivalent *linear* `ReactionSystem`, perform the Linear Mapping Approximation (LMA)
and return the corresponding linear [`RawMomentEquations`](@ref) of the system as well as
a Dictionary of reaction parameter substitutions obtained using LMA that are used to generate
the moment equations. See the [LMA theory section](@ref linear_mapping_approximation) for more details.

Notes:
- `rn_nonlinear` and `rn_linear` must be *identical* in layout in order to be interpreted correctly, and the
  nonlinear reactions contained in `rn_nonlinear` must all be linearised in `rn_linear` with rate coefficients
  updated accordingly. Although this requires a lot of manual input, automating the linearisation further is
  difficult due to arbitrary choices that may be mane in constructing the reaction networks.
- `binary_vars` *must* be specified for conditional closures as an array of indices of all species
  (as in [`speciesmap`](@ref)) which molecule number is a Bernoulli variable. Note that `rn_nonlinear`
  and `rn_linear` may internally order the species differently: `binary_vars` must be consistent with
  the ordering in the *nonlinear* network.
- By default the moment equations will be generated up to the order determined by the degree of nonlinearity
  of the nonlinear system's reactions. However, if higher order moment information is required, the optional
  `m_order` argument may be provided to increase the expansion order manually.
- `combinatoric_ratelaw=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
-  [`ReactionSystemMod`](@ref) is currently unsupported as by design it does not contain information of reaction
  substrates and hence requires an API change. To be updated in the future!
"""
function linear_mapping_approximation(rn_nonlinear::T, rn_linear::T, binary_vars::Array{Int,1}=Int[], m_order::Int=0;
                                      combinatoric_ratelaw = true) where T <: ReactionSystem

      # check that all necessary information is provided and that LMA is applicable on the given nonlinear system

      # TODO: add further checks that both rn_nonlinear and rn_linear have reactions ordered identically
      @assert !isempty(binary_vars) "LMA does not work if there are no binary species"
      submat = substoichmat(rn_nonlinear)'
      no_substrates = vcat(sum(submat, dims=2)...)
      no_binary = vcat(sum(submat[:, binary_vars], dims=2)...)
      nonlinear_rs_inds = findall(no_substrates .> 1)
      @assert all(x == 1 for x in no_binary[nonlinear_rs_inds]) "non-linear reactions must involve one type of binary species"
      @assert all(submat[nonlinear_rs_inds, binary_vars] .< 2) "cannot have more than 1 molecules of each binary species"
      @assert all(sum(substoichmat(rn_linear), dims=1) .<= 1) "the linear network cannot contain nonlinear reactions"

      linearised_rs = reactions(rn_linear)[nonlinear_rs_inds]
      @assert all(!reaction.only_use_rate for reaction in linearised_rs) "linearised nonlinear reactions must follow the law of mass action (defining reactions using Catalyst's →)"
      coeffs = (reaction.rate for reaction in linearised_rs)

      # apply LMA: substitute reaction parameters in the linear network moment equations to reflect
      # the nonlinear interactions according to the LMA methodology

      term_factors, term_powers, poly_order = polynomial_propensities(propensities(rn_nonlinear; combinatoric_ratelaw)[nonlinear_rs_inds], rn_nonlinear)

      order = iszero(m_order) ? poly_order : max(poly_order, m_order)
      try sys = generate_raw_moment_eqs(rn_linear, order; combinatoric_ratelaw, smap=speciesmap(rn_nonlinear))
      catch e
            error("LMA cannot handle reactions with non-polynomial rates\n $e")
      end
      sys = bernoulli_moment_eqs(sys, binary_vars)

      which_binary = [binary_vars[findfirst(submat[i, binary_vars] .> 0)] for i in nonlinear_rs_inds]
      sub_params = OrderedDict()
      μ = sys.μ
      for (coeff, factors, powers, binary_ind) in zip(coeffs, term_factors, term_powers, which_binary)
            sub_params[coeff] = sum(factor*μ[Tuple(power)] for (factor, power) in zip(factors, powers))
            #sub_params[coeff] /= μ[sys.iter_1[binary_ind]]
            sub_params[coeff] *= μ[sys.iter_1[binary_ind]]^-1
            sub_params[coeff] = simplify(sub_params[coeff])
      end

      LMA_eqs = Equation[]
      for eq in get_eqs(sys.odes)
            rhs = substitute(eq.rhs, sub_params)
            rhs = expand(rhs)
            push!(LMA_eqs, Equation(eq.lhs, rhs))
      end

      field_values = [getfield(sys, field) for field in fieldnames(typeof(sys))]

      iv = get_iv(sys.odes)
      ps = reactionparams(rn_nonlinear)
      vars = states(sys.odes)
      odename = Symbol(nameof(sys), "_LMA")
      odes = ODESystem(LMA_eqs, iv, vars, ps; name=odename)
      new_system = typeof(sys)(odes, field_values[2:end]...)

      new_system, sub_params

end

#= TODO: consider adding
  • steady-state stuff
  • identification of whether a linear system has a time-dependendent or steady-state solution
  • construction of probability distributions if a solution exists
=#
