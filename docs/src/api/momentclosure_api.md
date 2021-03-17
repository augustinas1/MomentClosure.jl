# [MomentClosure.jl API](@id api)
```@meta
CurrentModule = MomentClosure
```

## Model definition

MomentClosure is fully compatible with reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and stored as a [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Nevertheless, we have implemented a *limited* extension of `ModelingToolkit.ReactionSystem`, called `ReactionSystemMod`, that allows us to consider systems containing reactions whose products can be independent geometrically distributed random variables. For example, such reactions can be encountered in certain models of [autoregulatory gene networks](@ref gene_network_example). Moreover, `ReactionSystemMod` provides an alternative to Catalyst as the model is now constructed directly in terms of user-defined reaction propensity functions and the stoichiometric matrix, with system parameters and species defined using ModelingToolkit [`@parameters`](https://mtk.sciml.ai/stable/basics/ContextualVariables/#Parameters-1) and Symbolics [`@variables`](https://symbolics.juliasymbolics.org/stable/manual/variables/#Symbolics.@variables) respectively. However, we stress that defining models using Catalyst is preferred whenever possible due to its much richer integration within the broader [SciML](https://github.com/SciML/) ecosystem.

```@docs
@parameters
@variables
ReactionSystemMod
```

## Basic model properties

To keep the API consistent with Catalyst, we provide a number of identical [basic functions](https://catalyst.sciml.ai/stable/api/catalyst_api/#Basic-properties) allowing easy access to network properties:

```@docs
species
speciesmap
params
paramsmap
numspecies
numparams
numreactions
```

Moreover, we include two slight extensions of the Catalyst API:

```@docs
get_S_mat
propensities
```

## Moment Equations

```@docs
generate_raw_moment_eqs
RawMomentEquations
generate_central_moment_eqs
CentralMomentEquations
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
Note that the [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) (saved as `moment_eqs.odes`) can also be passed to `latexify` function directly.

Given [`ClosedMomentEquations`](@ref), the closure functions can be visualised in the same way by adding `:closure` argument:
```julia
latexify(moment_eqs, :closure)
```
## Stochastic Simulation Utilities

```@docs
get_raw_moments
get_central_moments
get_cumulants
```
