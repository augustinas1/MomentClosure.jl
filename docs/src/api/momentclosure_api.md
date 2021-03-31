# [MomentClosure.jl API](@id api)
```@meta
CurrentModule = MomentClosure
```

## Model definition

MomentClosure is fully compatible with reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and stored as a [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Nevertheless, we have implemented our own `ReactionSystemMod`, that allows us to consider systems containing reactions which products are independent geometrically distributed random variables (currently unsupported by Catalyst). For example, such reactions can be encountered in certain types of gene regulatory network models, as discussed in [this tutorial](@ref geometric-and-conditional). Moreover, `ReactionSystemMod` provides an alternative to Catalyst as the model can be now constructed directly in terms of user-defined reaction propensity functions and the stoichiometric matrix, with system parameters and species defined using ModelingToolkit [`@parameters`](https://mtk.sciml.ai/stable/basics/ContextualVariables/#Parameters-1) and Symbolics [`@variables`](https://symbolics.juliasymbolics.org/stable/manual/variables/#Symbolics.@variables) respectively.

We stress that [`ReactionSystemMod`](@ref), as the name suggests, is a rather trivial extension/modification of [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Although [`ReactionSystemMod`](@ref) provides the same basic functions to access network properties as does Catalyst for a `ModelingToolkit.ReactionSystem` ([see below](@ref api-basic-network-properties)), its integration within the broader SciML framework is very limited (only [compatible with SSA solvers](@ref stochastic_simulation_utilities)) and so we advise using Catalyst for model initialisation whenever possible. Finally, note that `ReactionSystemMod` is a temporary solution: we expect this interface to be supported by Catalyst directly soon (see the discussion [here](https://github.com/SciML/Catalyst.jl/issues/308) and [here](https://github.com/SciML/Catalyst.jl/issues/306)).

```@docs
@parameters
@variables
ReactionSystemMod
```

## [Basic model properties](@id api-basic-network-properties)

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

## [Stochastic Simulation Utilities](@id stochastic_simulation_utilities)

A [`ReactionSystemMod`](@ref) is compatible with the standard DifferentialEquations [`JumpProblem`](https://diffeq.sciml.ai/stable/types/jump_types/#Jump-Problems) and hence can be simulated using SSA. Internally, [`ReactionSystemMod`](@ref) is converted into a system of [ConstantRateJumps](https://diffeq.sciml.ai/stable/types/jump_types/#Types-of-Jumps:-Regular,-Variable,-Constant-Rate-and-Mass-Action) taking into account the geometrically distributed reaction products. Note, however, that the implementation here is very restricted: time-dependent propensities are currently not allowed and [MassActionJumps] (https://diffeq.sciml.ai/stable/types/jump_types/#Types-of-Jumps:-Regular,-Variable,-Constant-Rate-and-Mass-Action) cannot be used due to limitations of the [`ReactionSystemMod`](@ref) API, hence the performance is subpar.

We also provide provides functions for higher-order moment extraction from SSA data:
```@docs
get_raw_moments
get_central_moments
get_cumulants
```
