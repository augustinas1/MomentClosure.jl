# [MomentClosure.jl API](@id api)
```@meta
CurrentModule = MomentClosure
```

## Model definition

MomentClosure is fully compatible with reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and stored as a [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Nevertheless, we have implemented a *limited* extension of `ModelingToolkit.ReactionSystem`, called `ReactionSystemMod`, that allows us to consider systems containing reactions whose products can be independent geometrically distributed random variables. For example, such reactions can be encountered in certain models of [autoregulatory gene networks](gene_network_example.md). Moreover, `ReactionSystemMod` provides an alternative to Catalyst as the model is now constructed directly in terms of user-defined reaction propensity functions and the stoichiometric matrix, with system parameters and species defined using ModelingToolkit [`@parameters`](https://mtk.sciml.ai/stable/highlevel/#ModelingToolkit.@parameters) and [`@variables`](https://mtk.sciml.ai/stable/highlevel/#ModelingToolkit.@variables) respectively. However, we stress that defining models using Catalyst is preferred whenever possible due to its much richer integration within the broader [SciML](https://github.com/SciML/) ecosystem.

```@docs
ModelingToolkit.@parameters
ModelingToolkit.@variables
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


## Displaying equations and closures

```@docs
format_moment_eqs
format_closure
```
