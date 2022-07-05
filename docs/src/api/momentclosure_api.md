# [MomentClosure.jl API](@id api)
```@meta
CurrentModule = MomentClosure
```

## Model definition

MomentClosure is fully compatible with reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and stored as a [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Note that previously we had implemented our own `ReactionSystemMod`, that allowed us to consider systems containing reactions which products are independent geometrically distributed random variables. However, this is now deprecated as Catalyst has added support for such parameteric stoichiometries offering a much more complete and efficient feature set.

## [Basic model properties](@id api-basic-network-properties)

Moreover, we include a couple of tiny extensions to the Catalyst API:

```@docs
propensities
get_stoichiometry
```

## Moment Equations

```@docs
generate_raw_moment_eqs
RawMomentEquations
generate_central_moment_eqs
CentralMomentEquations
bernoulli_moment_eqs
```

## Moment Closure

```@docs
moment_closure
ClosedMomentEquations
deterministic_IC
```

## [Displaying Equations and Closures](@id visualisation_api)

The generated moment equations can be converted into LaTeX expressions using [Latexify](https://github.com/korsbo/Latexify.jl) as:
```julia
using Latexify
latexify(moment_eqs)
```
A [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) (saved as `moment_eqs.odes`) can also be passed to `latexify` function directly but the output will be different as we apply additional formatting to the symbolic expressions.

Given [`ClosedMomentEquations`](@ref), the closure functions can be visualised in the same way by adding a `:closure` argument:
```julia
latexify(moment_eqs, :closure)
```
Note that this will print out only those higher order moments which are found in the given moment equations. It is possible to print the closure functions of *all* higher order moments using `print_all=true` argument:
```julia
latexify(moment_eqs, :closure, print_all=true)
```

## Linear Mapping Approximation

```@docs
linear_mapping_approximation
```

## [Stochastic Simulation Utilities](@id stochastic_simulation_utilities)

We provide provides functions for higher-order moment extraction from SSA and FSP data:
```@docs
get_raw_moments
get_central_moments
get_cumulants
get_moments_FSP
```
