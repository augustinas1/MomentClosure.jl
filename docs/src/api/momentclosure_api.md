# MomentClosure.jl API
```@meta
CurrentModule = MomentClosure
```

## Model definition

As described in the [tutorial on model initialisation](@ref model_initialisation), MomentClosure is fully compatible with reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and stored as a [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Nevertheless, we have implemented a *limited* extension of `ModelingToolkit.ReactionSystem`, called `ReactionSystemMod`, that allows us to consider systems containing reactions whose products can be independent geometrically distributed random variables. Note that such systems, one being the gene network example [LINK], are currently not supported by Catalyst. Moreover, `ReactionSystemMod` provides an alternative to Catalyst as the model is now constructed directly in terms of user-defined reaction propensity functions and the stoichiometric matrix, with system parameters and species defined using ModelingToolkit [`@parameters`](https://mtk.sciml.ai/stable/highlevel/#ModelingToolkit.@parameters) and [`@variables`](https://mtk.sciml.ai/stable/highlevel/#ModelingToolkit.@variables) respectively. However, we stress that defining models using Catalyst is preferred whenever possible due to its much richer integration within the broader [SciML](https://github.com/SciML/) ecosystem.

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
