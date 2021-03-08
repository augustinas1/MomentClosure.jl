# [Using MomentClosure](@id main_tutorial)

This tutorial is an introduction to using MomentClosure to define chemical reaction network models, generate the corresponding moment equations, apply moment closure approximations and finally solve the resulting system of ODEs. To demonstrate this functionality, we will consider a specific case of an oscillatory chemical system known as the Brusselator, characterised by the reactions [1]
```math
\begin{align*}
2X + Y &\stackrel{c_1}{\rightarrow} 3X, \\
X &\stackrel{c_2}{\rightarrow} Y, \\
∅ &\underset{c_4}{\stackrel{c_3}{\rightleftharpoons}} X.
\end{align*}
```
We have chosen this particular model to start with as it has been studied with different moment closures before by Schnoerr et al. (2015) [1] and so it is useful as a reference point. The plots of moment trajectories we obtain in this tutorial fully reproduce some of the figures published in the paper, hence (partially) proving the validity of our implementation.

The terminology and the notation used throughout is consistent with the **Theory** section of the docs and we advise giving it a skim-through.

## Model Initialisation

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) provides a comprehensive interface to modelling chemical reaction networks in Julia and can be used to construct models fully-compatible with MomentClosure. For more details on how to do so we recommend reading [Catalyst's tutorial](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/). This way, the Brusselator can be defined as:
```julia
using Catalyst
rn = @reaction_network begin
  (c₁/Ω^2), 2X + Y → 3X
  (c₂), X → Y
  (c₃*Ω, c₄), 0 ↔ X
end c₁ c₂ c₃ c₄ Ω
```
The returned `rn` is an instance of [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). Note that we have also included the system-size parameter $\Omega$ which will be of interest later on looking at the system's dynamics.

**Alternatively**, models can be defined as a MomentClosure's built-in [`ReactionSystemMod`](@ref) by considering the net stoichiometry matrix, $S$, and a vector of the corresponding reaction propensity functions, $\mathbf{a}$. As in Catalyst, the model specification is based on the rich symbolic-numeric modelling framework provided by [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). Note that we have to explicitly define the molecule numbers of each chemical species as [`Symbolics.@variables`](@ref) and reaction constants as [`ModelingToolkit.@parameters`](@ref). The time, $t$, must also be initialised as a `ModelingToolkit.parameter` and used to indicate the time-dependence of species' variables. The reaction network can then be constructed as follows:
```julia
using MomentClosure
@parameters t, c₁, c₂, c₃, c₄, Ω
@variables X(t), Y(t)

# stoichiometric matrix
S_mat = [ 1 -1  1 -1;
         -1  1  0  0]

# propensity functions
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]

rn2 = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)
```
We stress that [`ReactionSystemMod`](@ref), as the name suggests, is a rather trivial extension/modification of [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem). The main motivation behind it is to add support for systems that include reactions whose products are independent geometrically distributed random variables ([currently unsupported by Catalyst](https://github.com/SciML/Catalyst.jl/issues/308)), see [the auto-regulatory gene network example](gene_network_example.md) for more details. Although [`ReactionSystemMod`](@ref) provides the same basic functions to access network properties as does Catalyst for a `ModelingToolkit.ReactionSystem` ([see the API](@ref api)), its functionality is nowhere near as rich. For example, `ReactionSystemMod` cannot be used directly for [deterministic or stochastic simulations](https://github.com/SciML/DifferentialEquations.jl/) using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl/). Therefore, we advise using Catalyst for model initialisation whenever possible.

Note that the net stoichiometry matrix and an array of the corresponding propensities, if needed, can be extracted directly from a `ModelingToolkit.ReactionSystem` using functions [`get_S_mat`](@ref) and [`propensities`](@ref) respectively.

## Generating Moment Equations

Having defined the Brusselator model, we can proceed to obtain the moment equations. The system follows the law of mass action, i.e., all propensity functions are polynomials in molecule numbers $X(t)$ and $Y(t)$, and so we can generate either *raw* or *central* moment equations, as described in the **Theory** section on [moment expansion](@ref moment_expansion).

#### Raw Moment Equations

Let's start with the raw moment equations which we choose to generate up to second order ($m=2$):
```julia
using MomentClosure
raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=false)
```
Note that we have set `combinatoric_ratelaw=false` in order to ignore the factorial scaling factors which [Catalyst adds to mass-action reactions](https://catalyst.sciml.ai/stable/tutorials/models/#Reaction-rate-laws-used-in-simulations). The function [`generate_raw_moment_eqs`](@ref) returns an instance of [`RawMomentEquations`](@ref MomentClosure.RawMomentEquations) that contains a [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) composed of all the moment equations (accessed by `raw_eqs.odes`).

We can use [Latexify](https://github.com/korsbo/Latexify.jl) to look at the generated moment equations:
```julia
using Latexify
latexify(raw_eqs)
```
```math
\begin{align*}
\frac{d\mu{_{10}}}{dt} =& c{_3} \Omega + c{_1} \mu{_{21}} \Omega^{-2} - c{_2} \mu{_{10}} - c{_4} \mu{_{10}} - c{_1} \mu{_{11}} \Omega^{-2} \\
\frac{d\mu{_{01}}}{dt} =& c{_2} \mu{_{10}} + c{_1} \mu{_{11}} \Omega^{-2} - c{_1} \mu{_{21}} \Omega^{-2} \\
\frac{d\mu{_{20}}}{dt} =& c{_2} \mu{_{10}} + c{_3} \Omega + c{_4} \mu{_{10}} + 2 c{_1} \mu{_{31}} \Omega^{-2} + 2 c{_3} \Omega \mu{_{10}} - 2 c{_2} \mu{_{20}} - 2 c{_4} \mu{_{20}} - c{_1} \mu{_{11}} \Omega^{-2} - c{_1} \mu{_{21}} \Omega^{-2} \\
\frac{d\mu{_{11}}}{dt} =& c{_2} \mu{_{20}} + c{_1} \mu{_{11}} \Omega^{-2} + c{_1} \mu{_{22}} \Omega^{-2} + c{_3} \Omega \mu{_{01}} - c{_2} \mu{_{10}} - c{_2} \mu{_{11}} - c{_4} \mu{_{11}} - c{_1} \mu{_{12}} \Omega^{-2} - c{_1} \mu{_{31}} \Omega^{-2} \\
\frac{d\mu{_{02}}}{dt} =& c{_2} \mu{_{10}} + c{_1} \mu{_{21}} \Omega^{-2} + 2 c{_2} \mu{_{11}} + 2 c{_1} \mu{_{12}} \Omega^{-2} - c{_1} \mu{_{11}} \Omega^{-2} - 2 c{_1} \mu{_{22}} \Omega^{-2}
\end{align*}
```
The raw moments are defined as
```math
\mu_{ij}(t) = \langle X(t)^i Y(t)^j \rangle
```
where $\langle \rangle$ denote the expectation value and we have explicitly included the time-dependence for completeness (made implicit in the formatted moment equations). Note that the ordering of species ($X$ first and $Y$ second) is consistent with the order these variables appear within the [`Catalyst.@reaction_network`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.@reaction_network) macro. The ordering can also be checked using [`Catalyst.speciesmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.speciesmap) function:
```julia
speciesmap(rn)
```
```julia
Dict{Term{Real},Int64} with 2 entries:
  X(t) => 1
  Y(t) => 2
```
Note that [`speciesmap`](@ref) can be used in the same way with a [`ReactionSystemMod`](@ref).

Coming back to the generated moment equations, we observe that they depend on higher-order moments. For example, the ODE for $\mu_{02}$ depends on third order moments $μ_{12}$ and $μ_{21}$ and the fourth order moment $\mu_{22}$. Consider the general case of [raw moment equations](@ref raw_moment_eqs): if a network involves reactions that are polynomials (in molecule numbers) of *at most* order $k$, then its $m^{\text{th}}$ order moment equations will depend on moments up to order $m+k-1$. Hence the relationship seen above is expected as the Brusselator involves a trimolecular reaction whose corresponding propensity function is a third order polynomial in $X(t)$ and $Y(t)$. The number denoting the highest order of moments encountered in the generated [`RawMomentEquations`](@ref MomentClosure.RawMomentEquations) can also be accessed as `raw_eqs.q_order` (returning `4` in this case).

#### Central Moment Equations

The corresponding central moment equations can also be easily generated:
```julia
central_eqs = generate_central_moment_eqs(rn, 2, combinatoric_ratelaw=false)
```
Note that in case of non-polynomial propensity functions the [expansion order $q$](@ref central_moment_eqs) must also be specified, see the [P53 system example](P53_system_example.md) for more details. Luckily, the Brusselator contains only mass-action reactions and hence $q$ is automatically determined by the highest order (polynomial) propensity. The function [`generate_central_moment_eqs`](@ref) returns an instance of [`CentralMomentEquations`](@ref MomentClosure.CentralMomentEquations). As before, we can visualise the central moment equations:
```julia
latexify(central_eqs)
```
```math
\begin{align*}
\frac{d\mu{_{10}}}{dt} =& c{_3} \Omega + c{_1} M{_{21}} \Omega^{-2} + c{_1} M{_{20}} \mu{_{01}} \Omega^{-2} + c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} + 2 c{_1} M{_{11}} \mu{_{10}} \Omega^{-2} - c{_2} \mu{_{10}} - c{_4} \mu{_{10}} - c{_1} M{_{11}} \Omega^{-2} - c{_1} \mu{_{01}} \mu{_{10}} \Omega^{-2} \\
\frac{d\mu{_{01}}}{dt} =& c{_2} \mu{_{10}} + c{_1} M{_{11}} \Omega^{-2} + c{_1} \mu{_{01}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{21}} \Omega^{-2} - 2 c{_1} M{_{11}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{20}} \mu{_{01}} \Omega^{-2} - c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} \\
\frac{dM{_{20}}}{dt} =& c{_2} \mu{_{10}} + c{_3} \Omega + c{_4} \mu{_{10}} + 2 c{_1} M{_{31}} \Omega^{-2} + c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} + 2 c{_1} M{_{11}} \Omega^{-2} \mu{_{10}}^{2} + 4 c{_1} M{_{21}} \mu{_{10}} \Omega^{-2} + 2 c{_1} M{_{30}} \mu{_{01}} \Omega^{-2} + 4 c{_1} M{_{20}} \mu{_{01}} \mu{_{10}} \Omega^{-2} - 2 c{_2} M{_{20}} - 2 c{_4} M{_{20}} - c{_1} M{_{11}} \Omega^{-2} - c{_1} M{_{21}} \Omega^{-2} - c{_1} M{_{20}} \mu{_{01}} \Omega^{-2} - c{_1} \mu{_{01}} \mu{_{10}} \Omega^{-2} \\
\frac{dM{_{11}}}{dt} =& c{_2} M{_{20}} + c{_1} M{_{11}} \Omega^{-2} + c{_1} M{_{22}} \Omega^{-2} + c{_1} M{_{02}} \Omega^{-2} \mu{_{10}}^{2} + c{_1} M{_{21}} \mu{_{01}} \Omega^{-2} + c{_1} \mu{_{01}} \mu{_{10}} \Omega^{-2} + 2 c{_1} M{_{12}} \mu{_{10}} \Omega^{-2} + 2 c{_1} M{_{11}} \mu{_{01}} \mu{_{10}} \Omega^{-2} - c{_2} M{_{11}} - c{_2} \mu{_{10}} - c{_4} M{_{11}} - c{_1} M{_{12}} \Omega^{-2} - c{_1} M{_{31}} \Omega^{-2} - c{_1} M{_{02}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{11}} \mu{_{01}} \Omega^{-2} - c{_1} M{_{11}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{11}} \Omega^{-2} \mu{_{10}}^{2} - 2 c{_1} M{_{21}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{30}} \mu{_{01}} \Omega^{-2} - c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} - 2 c{_1} M{_{20}} \mu{_{01}} \mu{_{10}} \Omega^{-2} \\
\frac{dM{_{02}}}{dt} =& c{_2} \mu{_{10}} + c{_1} M{_{21}} \Omega^{-2} + 2 c{_2} M{_{11}} + 2 c{_1} M{_{12}} \Omega^{-2} + c{_1} M{_{20}} \mu{_{01}} \Omega^{-2} + c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} + 2 c{_1} M{_{02}} \mu{_{10}} \Omega^{-2} + 2 c{_1} M{_{11}} \mu{_{01}} \Omega^{-2} + 2 c{_1} M{_{11}} \mu{_{10}} \Omega^{-2} - c{_1} M{_{11}} \Omega^{-2} - 2 c{_1} M{_{22}} \Omega^{-2} - 2 c{_1} M{_{02}} \Omega^{-2} \mu{_{10}}^{2} - 4 c{_1} M{_{12}} \mu{_{10}} \Omega^{-2} - 2 c{_1} M{_{21}} \mu{_{01}} \Omega^{-2} - c{_1} \mu{_{01}} \mu{_{10}} \Omega^{-2} - 4 c{_1} M{_{11}} \mu{_{01}} \mu{_{10}} \Omega^{-2}
\end{align*}
```
Unfortunately, central moment equations often take a visually painful form. Note that the first two ODEs, as before, indicate the means, and the central moments are denoted as
```math
M_{ij}(t) = \langle (X(t)-\mu_{10}(t))^i (Y(t)-\mu_{01}(t))^j \rangle.
```

## Applying Moment Closure

As observed above, the moment equations of the Brusselator are coupled and depend on higher order moments—we have an infinite hierarchy of ODEs in our hands which cannot be solved directly and requires approximate treatment. One way of approaching the problem is to apply moment closure approximations (MAs), in which higher order moments are expressed as functions of lower order moments, thus effectively truncating the hierarchy and enabling a numerical solution. A variety of MAs have been proposed in literature and are implemented in MomentClosure.jl, see the **Theory** [section on MAs](@ref moment_closure_approximations) for more details.

Let's apply [normal closure](@ref normal_closure) to the raw moment equations `raw_eqs` we have generated earlier using [`moment_closure`](@ref) function:
```julia
closed_raw_eqs = moment_closure(raw_eqs, "normal")
```
The function returns [`ClosedMomentEquations`](@ref MomentClosure.ClosedMomentEquations) that consists of both the [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) containing all closed moment equations as well as the specific closure functions for each higher order raw moments. We can use Latexify again to look at the closed ODEs:
```julia
latexify(closed_raw_eqs)
```
```math
\begin{align*}
\frac{d\mu{_{10}}}{dt} =& c{_3} \Omega + c{_1} \mu{_{01}} \mu{_{20}} \Omega^{-2} + 2 c{_1} \mu{_{10}} \mu{_{11}} \Omega^{-2} - c{_2} \mu{_{10}} - c{_4} \mu{_{10}} - c{_1} \mu{_{11}} \Omega^{-2} - 2 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} \\
\frac{d\mu{_{01}}}{dt} =& c{_2} \mu{_{10}} + c{_1} \mu{_{11}} \Omega^{-2} + 2 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} - c{_1} \mu{_{01}} \mu{_{20}} \Omega^{-2} - 2 c{_1} \mu{_{10}} \mu{_{11}} \Omega^{-2} \\
\frac{d\mu{_{20}}}{dt} =& c{_2} \mu{_{10}} + c{_3} \Omega + c{_4} \mu{_{10}} + 2 c{_3} \Omega \mu{_{10}} + 2 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} + 6 c{_1} \mu{_{11}} \mu{_{20}} \Omega^{-2} - 2 c{_2} \mu{_{20}} - 2 c{_4} \mu{_{20}} - c{_1} \mu{_{11}} \Omega^{-2} - c{_1} \mu{_{01}} \mu{_{20}} \Omega^{-2} - 4 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{3} - 2 c{_1} \mu{_{10}} \mu{_{11}} \Omega^{-2} \\
\frac{d\mu{_{11}}}{dt} =& c{_2} \mu{_{20}} + c{_1} \mu{_{11}} \Omega^{-2} + c{_3} \Omega \mu{_{01}} + c{_1} \mu{_{02}} \mu{_{20}} \Omega^{-2} + 2 c{_1} \Omega^{-2} \mu{_{11}}^{2} + 2 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{3} + 2 c{_1} \mu{_{10}} \Omega^{-2} \mu{_{01}}^{2} - c{_2} \mu{_{10}} - c{_2} \mu{_{11}} - c{_4} \mu{_{11}} - 2 c{_1} \mu{_{01}} \mu{_{11}} \Omega^{-2} - c{_1} \mu{_{02}} \mu{_{10}} \Omega^{-2} - 3 c{_1} \mu{_{11}} \mu{_{20}} \Omega^{-2} - 2 c{_1} \Omega^{-2} \mu{_{01}}^{2} \mu{_{10}}^{2} \\
\frac{d\mu{_{02}}}{dt} =& c{_2} \mu{_{10}} + 2 c{_2} \mu{_{11}} + c{_1} \mu{_{01}} \mu{_{20}} \Omega^{-2} + 4 c{_1} \mu{_{01}} \mu{_{11}} \Omega^{-2} + 2 c{_1} \mu{_{02}} \mu{_{10}} \Omega^{-2} + 2 c{_1} \mu{_{10}} \mu{_{11}} \Omega^{-2} + 4 c{_1} \Omega^{-2} \mu{_{01}}^{2} \mu{_{10}}^{2} - c{_1} \mu{_{11}} \Omega^{-2} - 4 c{_1} \Omega^{-2} \mu{_{11}}^{2} - 2 c{_1} \mu{_{01}} \Omega^{-2} \mu{_{10}}^{2} - 2 c{_1} \mu{_{02}} \mu{_{20}} \Omega^{-2} - 4 c{_1} \mu{_{10}} \Omega^{-2} \mu{_{01}}^{2}
\end{align*}
```
The closure functions can also be displayed by adding `:closure` argument:
```julia
latexify(exprs, :closure)
```
```math
\begin{align*}
\mu{_{30}} =& 3 \mu{_{10}} \mu{_{20}} - 2 \mu{_{10}}^{3} \\
\mu{_{21}} =& \mu{_{01}} \mu{_{20}} + 2 \mu{_{10}} \mu{_{11}} - 2 \mu{_{01}} \mu{_{10}}^{2} \\
\mu{_{12}} =& \mu{_{02}} \mu{_{10}} + 2 \mu{_{01}} \mu{_{11}} - 2 \mu{_{10}} \mu{_{01}}^{2} \\
\mu{_{03}} =& 3 \mu{_{01}} \mu{_{02}} - 2 \mu{_{01}}^{3} \\
\mu{_{40}} =& 6 \mu{_{10}}^{4} + 3 \mu{_{20}}^{2} + 4 \mu{_{10}} \mu{_{30}} - 12 \mu{_{20}} \mu{_{10}}^{2} \\
\mu{_{31}} =& \mu{_{01}} \mu{_{30}} + 6 \mu{_{01}} \mu{_{10}}^{3} + 3 \mu{_{10}} \mu{_{21}} + 3 \mu{_{11}} \mu{_{20}} - 6 \mu{_{11}} \mu{_{10}}^{2} - 6 \mu{_{01}} \mu{_{10}} \mu{_{20}} \\
\mu{_{22}} =& \mu{_{02}} \mu{_{20}} + 2 \mu{_{11}}^{2} + 2 \mu{_{01}} \mu{_{21}} + 2 \mu{_{10}} \mu{_{12}} + 6 \mu{_{01}}^{2} \mu{_{10}}^{2} - 2 \mu{_{02}} \mu{_{10}}^{2} - 2 \mu{_{20}} \mu{_{01}}^{2} - 8 \mu{_{01}} \mu{_{10}} \mu{_{11}} \\
\mu{_{13}} =& \mu{_{03}} \mu{_{10}} + 3 \mu{_{01}} \mu{_{12}} + 3 \mu{_{02}} \mu{_{11}} + 6 \mu{_{10}} \mu{_{01}}^{3} - 6 \mu{_{11}} \mu{_{01}}^{2} - 6 \mu{_{01}} \mu{_{02}} \mu{_{10}} \\
\mu{_{04}} =& 6 \mu{_{01}}^{4} + 3 \mu{_{02}}^{2} + 4 \mu{_{01}} \mu{_{03}} - 12 \mu{_{02}} \mu{_{01}}^{2}
\end{align*}
```
Similarly, we can close central moment equations using
```julia
closed_central_eqs = moment_closure(central_eqs, "normal")
```
and print out the corresponding closure functions:
```julia
latexify(closed_central_eqs, :closure)
```
```math
\begin{align*}
M{_{30}} =& 0 \\
M{_{21}} =& 0 \\
M{_{12}} =& 0 \\
M{_{03}} =& 0 \\
M{_{40}} =& 3 M{_{20}}^{2} \\
M{_{31}} =& 3 M{_{11}} M{_{20}} \\
M{_{22}} =& M{_{02}} M{_{20}} + 2 M{_{11}}^{2} \\
M{_{13}} =& 3 M{_{02}} M{_{11}} \\
M{_{04}} =& 3 M{_{02}}^{2}
\end{align*}
```
Higher order central moments under normal closure take a rather simple form compared to their raw moment equivalents. That can be expected due to the relationship between central moments and cumulants, on which [the closure is based](@ref normal_closure).

## Solving Moment Equations

The closed moment equations can be solved numerically using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl/) that provides a variety of highly-efficient solvers and analysis tools. In order to do so, first we need to specify the values of all system parameters, the initial condition and the timespan to solve over. Then the [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) corresponding to the moment equations can be [directly converted](https://mtk.sciml.ai/stable/#Compatible-Numerical-Solvers-1) into an `ODEProblem` which can finally be solved. Let's go through the procedure step-by-step for the closed raw moment equations (`closed_raw_eqs`).

Note that

First we need to create a map between the symbolic parameter variables and their specified values.
Knowing the ordering of parameters (accessible using `paramsmap` function) we can do that as follows:
```julia
# parameters [c₁, c₂, c₃, c₄, Ω]
p = [0.9, 2, 1, 1, 100]
pmap = Pair.(params(rn), p)
```
[See the Catalyst tutorial for a closely related example](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/#Mass-Action-ODE-Models). As an alternative, we can construct `pmap` manually:
```julia
@parameters c₁, c₂, c₃, c₄, Ω # if not in scope
pmap = [c₁ => 0.9,
        c₂ => 2,
        c₃ => 1,
        c₄ => 1,
        Ω => 100]
```
Next we can specify the initial condition. Usually with moment equations we are considering *deterministic* initial conditions so that the molecule numbers at initial time take the specified integer values with probability one. Hence in our case we set $X(t=0) = X_0$ and $Y(t=0) = Y_0$. Probability one implies that initially the means will be equal to the molecule numbers, i.e., $μ_{10}(t=0) = X_0$ and $μ_{01}(t=0) = Y_0$, whereas all higher order raw moments will be products of the corresponding powers of the means, e.g., $μ_{21} = X_0^2 Y_0$. Note that all central moments would be set to zero in this case. To make life easier we use [`deterministic_IC`](@ref) function which, given the initial molecule numbers, automatically constructs the variable mapping under deterministic initial conditions:
```julia
# initial molecule numbers [X, Y] (order as in speciesmap(rn))
u₀ = [0, 0]
u₀map = deterministic_IC(u₀, closed_raw_eqs)
```
The next ingredient, the time interval to solve on, can be specified simply as:
```julia
tspan = (0., 100.)
```
Now we can create the corresponding `ODEProblem`:
```julia
using OrdinaryDiffEq
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, pmap)
```
Note that we are using only [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) as [there is no need to load the whole DifferentialEquations library](https://diffeq.sciml.ai/stable/features/low_dep/#Low-Dependency-Usage).

Finally, we have everything we need to solve the raw moment equations which can be done using any ODE solver [implemented within DifferentialEquations.jl](https://diffeq.sciml.ai/dev/solvers/ode_solve/). We use the default `Tsit5()` solver and then [plot](https://diffeq.sciml.ai/stable/basics/plot/#plot) the obtained mean molecule numbers:
```julia
sol = solve(oprob, Tsit5(), saveat=0.1)

using Plots
plot(sol, vars=(0, [1,2]), lw=2)
```
![Brusselator means 1](../assets/brusselator_means_1.svg)

The obtained moment dynamics show damped oscillations which is the expected averaged behaviour of the Brusselator in a vast swathe of parameter space [1]. However, to establish more clearly how well a second order moment expansion with normal moment closure performs for this system and this specific set of parameters, we can compare the resulting moment trajectories to the true moment estimates obtained using [Gillespie's Stochastic Simulation Algorithm (SSA)](https://en.wikipedia.org/wiki/Gillespie_algorithm).

## Stochastic Simulation


## Common Problems

## References

[1]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Comparison of different moment-closure approximations for stochastic chemical kinetics", The Journal of Chemical Physics 143, 185101 (2015). [https://doi.org/10.1063/1.4934990](https://doi.org/10.1063/1.4934990)
