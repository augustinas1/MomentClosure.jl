# [Using MomentClosure](@id main_tutorial)

This tutorial is an introduction to using MomentClosure to define chemical reaction network models, generate the corresponding moment equations, apply moment closure approximations and finally solve the resulting system of ODEs. To demonstrate this functionality, we will consider a specific case of an oscillatory chemical system known as the Brusselator, characterised by the reactions
```math
\begin{align*}
2X + Y &\stackrel{c_1}{\rightarrow} 3X, \\
X &\stackrel{c_2}{\rightarrow} Y, \\
∅ &\underset{c_4}{\stackrel{c_3}{\rightleftharpoons}} X.
\end{align*}
```
We have chosen this particular model to start with as it has been studied with different moment closures before by [Schnoerr et al. (2015)](https://doi.org/10.1063/1.4934990) and so it is useful as a reference point. The plots of moment trajectories we obtain in this tutorial fully reproduce some of the figures published in the paper, hence (partially) proving the validity of our implementation.

The terminology and the notation used throughout is consistent with the **Theory** section of the docs and we advise giving it a skim-through.

## Model Initialisation

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) provides a comprehensive interface to modelling chemical reaction networks in Julia and can be used to construct models fully compatible with MomentClosure. For more details on how to do so we recommend reading [Catalyst's tutorial](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/introduction_to_catalyst/). This way, the Brusselator can be defined as:
```julia
using Catalyst
rn = @reaction_network begin
  # including system-size parameter Ω
  @parameters c₁ c₂ c₃ c₄ Ω
  (c₁/Ω^2), 2X + Y → 3X
  (c₂), X → Y
  (c₃*Ω, c₄), 0 ↔ X
end
```
The returned `rn` is an instance of [`Catalyst.ReactionSystem`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.ReactionSystem). The net stoichiometry matrix and an array of the corresponding propensities, if needed, can be extracted directly from the model using [`Catalyst.netstoichmat`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.netstoichmat) and MomentClosure function [`propensities`](@ref) respectively.

Note that MomentClosure also supports systems containing geometrically distributed reaction products that can be defined [using Catalyst](https://docs.sciml.ai/Catalyst/stable/model_creation/parametric_stoichiometry/)—see [this tutorial](@ref geometric-and-conditional) for more details.

## Generating Moment Equations

We can now obtain the moment equations. The system follows the law of mass action, i.e., all propensity functions are polynomials in molecule numbers $X(t)$ and $Y(t)$, and so we can generate either *raw* or *central* moment equations, as described in the **Theory** section on [moment expansion](@ref moment_expansion_CME).

#### Raw Moment Equations

Let's start with the raw moment equations which we choose to generate up to second order ($m=2$):
```julia
using MomentClosure
raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaws=false)
```
Note that we have set `combinatoric_ratelaws=false` in order to ignore the factorial scaling factors which [Catalyst adds to mass-action reactions](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/introduction_to_catalyst/#introduction_to_catalyst_ratelaws). The function [`generate_raw_moment_eqs`](@ref) returns an instance of [`RawMomentEquations`](@ref MomentClosure.RawMomentEquations) that contains a [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) composed of all the moment equations (accessed by `raw_eqs.odes`).

We can use [Latexify](https://github.com/korsbo/Latexify.jl) to look at the generated moment equations:
```julia
using Latexify
latexify(raw_eqs)
```
```math
\begin{align*}
\frac{d\mu_{1 0}}{dt} =& c_3 \Omega + c_1 \Omega^{-2} \mu_{2 1} - c_2 \mu_{1 0} - c_4 \mu_{1 0} - c_1 \Omega^{-2} \mu_{1 1} \\
\frac{d\mu_{0 1}}{dt} =& c_2 \mu_{1 0} + c_1 \Omega^{-2} \mu_{1 1} - c_1 \Omega^{-2} \mu_{2 1} \\
\frac{d\mu_{2 0}}{dt} =& c_2 \mu_{1 0} + c_3 \Omega + c_4 \mu_{1 0} + 2 c_3 \Omega \mu_{1 0} + 2 c_1 \Omega^{-2} \mu_{3 1} - 2 c_2 \mu_{2 0} - 2 c_4 \mu_{2 0} - c_1 \Omega^{-2} \mu_{1 1} - c_1 \Omega^{-2} \mu_{2 1} \\
\frac{d\mu_{1 1}}{dt} =& c_2 \mu_{2 0} + c_3 \Omega \mu_{0 1} + c_1 \Omega^{-2} \mu_{1 1} + c_1 \Omega^{-2} \mu_{2 2} - c_2 \mu_{1 0} - c_2 \mu_{1 1} - c_4 \mu_{1 1} - c_1 \Omega^{-2} \mu_{1 2} - c_1 \Omega^{-2} \mu_{3 1} \\
\frac{d\mu_{0 2}}{dt} =& c_2 \mu_{1 0} + 2 c_2 \mu_{1 1} + c_1 \Omega^{-2} \mu_{2 1} + 2 c_1 \Omega^{-2} \mu_{1 2} - c_1 \Omega^{-2} \mu_{1 1} - 2 c_1 \Omega^{-2} \mu_{2 2}
\end{align*}
```
The raw moments are defined as
```math
\mu_{ij}(t) = \langle X(t)^i Y(t)^j \rangle
```
where $\langle \rangle$ denote the expectation value and we have explicitly included the time-dependence for completeness (made implicit in the formatted moment equations). Note that the ordering of species ($X$ first and $Y$ second) is consistent with the order these variables appear within the [`Catalyst.@reaction_network`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.@reaction_network) macro. The ordering can also be checked using [`Catalyst.speciesmap`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.speciesmap) function:
```julia
speciesmap(raw_eqs)
```
```julia
Dict{Term{Real},Int64} with 2 entries:
  X(t) => 1
  Y(t) => 2
```

Coming back to the generated moment equations, we observe that they depend on higher-order moments. For example, the ODE for $\mu_{02}$ depends on third order moments $μ_{12}$ and $μ_{21}$ and the fourth order moment $\mu_{22}$. Consider the general case of [raw moment equations](@ref raw_moment_eqs): if a network involves reactions that are polynomials (in molecule numbers) of *at most* order $k$, then its $m^{\text{th}}$ order moment equations will depend on moments up to order $m+k-1$. Hence the relationship seen above is expected as the Brusselator involves a trimolecular reaction whose corresponding propensity function is a third order polynomial in $X(t)$ and $Y(t)$. The number denoting the highest order of moments encountered in the generated [`RawMomentEquations`](@ref MomentClosure.RawMomentEquations) can also be accessed as `raw_eqs.q_order` (returning `4` in this case).

#### Central Moment Equations

The corresponding central moment equations can also be easily generated:
```julia
central_eqs = generate_central_moment_eqs(rn, 2, combinatoric_ratelaws=false)
```
Note that in case of non-polynomial propensity functions the [Taylor expansion order $q$](@ref central_moment_eqs) must also be specified, see the [P53 system example](P53_system_example.md) for more details. Luckily, the Brusselator contains only mass-action reactions and hence $q$ is automatically determined by the highest order (polynomial) propensity. The function [`generate_central_moment_eqs`](@ref) returns an instance of [`CentralMomentEquations`](@ref MomentClosure.CentralMomentEquations). As before, we can visualise the central moment equations:
```julia
latexify(central_eqs)
```
```math
\begin{align*}
\frac{d\mu_{1 0}}{dt} =& c_3 \Omega + c_1 \Omega^{-2} M_{2 1} + c_1 \Omega^{-2} \mu_{1 0}^{2} \mu_{0 1} + c_1 \Omega^{-2} M_{2 0} \mu_{0 1} + 2 c_1 \Omega^{-2} M_{1 1} \mu_{1 0} - c_2 \mu_{1 0} - c_4 \mu_{1 0} - c_1 \Omega^{-2} M_{1 1} - c_1 \Omega^{-2} \mu_{0 1} \mu_{1 0} \\
\frac{d\mu_{0 1}}{dt} =& c_2 \mu_{1 0} + c_1 \Omega^{-2} M_{1 1} + c_1 \Omega^{-2} \mu_{0 1} \mu_{1 0} - c_1 \Omega^{-2} M_{2 1} - c_1 \Omega^{-2} \mu_{1 0}^{2} \mu_{0 1} - c_1 \Omega^{-2} M_{2 0} \mu_{0 1} - \frac{2}{1} c_1 \Omega^{-2} M_{1 1} \mu_{1 0} \\
\frac{dM_{2 0}}{dt} =& c_3 \Omega + c_2 \mu_{1 0} + c_4 \mu_{1 0} + 2 c_1 \Omega^{-2} M_{3 1} + c_1 \Omega^{-2} \mu_{1 0}^{2} \mu_{0 1} + 2 c_1 \Omega^{-2} \mu_{1 0}^{2} M_{1 1} + 2 c_1 \Omega^{-2} M_{3 0} \mu_{0 1} + 4 c_1 \Omega^{-2} M_{2 1} \mu_{1 0} + 4 c_1 \Omega^{-2} M_{2 0} \mu_{0 1} \mu_{1 0} - 2 c_2 M_{2 0} - 2 c_4 M_{2 0} - c_1 \Omega^{-2} M_{1 1} - c_1 \Omega^{-2} M_{2 1} - c_1 \Omega^{-2} M_{2 0} \mu_{0 1} - c_1 \Omega^{-2} \mu_{0 1} \mu_{1 0} \\
\frac{dM_{1 1}}{dt} =& c_2 M_{2 0} + c_1 \Omega^{-2} M_{1 1} + c_1 \Omega^{-2} M_{2 2} + c_1 \Omega^{-2} \mu_{1 0}^{2} M_{0 2} + c_1 \Omega^{-2} M_{2 1} \mu_{0 1} + c_1 \Omega^{-2} \mu_{0 1} \mu_{1 0} + 2 c_1 \Omega^{-2} M_{1 2} \mu_{1 0} + 2 c_1 \Omega^{-2} M_{1 1} \mu_{0 1} \mu_{1 0} - c_2 M_{1 1} - c_4 M_{1 1} - c_2 \mu_{1 0} - c_1 \Omega^{-2} M_{1 2} - c_1 \Omega^{-2} M_{3 1} - c_1 \Omega^{-2} \mu_{1 0}^{2} M_{1 1} - c_1 \Omega^{-2} \mu_{1 0}^{2} \mu_{0 1} - c_1 \Omega^{-2} M_{0 2} \mu_{1 0} - c_1 \Omega^{-2} M_{1 1} \mu_{0 1} - c_1 \Omega^{-2} M_{3 0} \mu_{0 1} - c_1 \Omega^{-2} M_{1 1} \mu_{1 0} - 2 c_1 \Omega^{-2} M_{2 1} \mu_{1 0} - 2 c_1 \Omega^{-2} M_{2 0} \mu_{0 1} \mu_{1 0} \\
\frac{dM_{0 2}}{dt} =& c_2 \mu_{1 0} + 2 c_2 M_{1 1} + c_1 \Omega^{-2} M_{2 1} + 2 c_1 \Omega^{-2} M_{1 2} + c_1 \Omega^{-2} \mu_{1 0}^{2} \mu_{0 1} + c_1 \Omega^{-2} M_{2 0} \mu_{0 1} + 2 c_1 \Omega^{-2} M_{0 2} \mu_{1 0} + 2 c_1 \Omega^{-2} M_{1 1} \mu_{0 1} + 2 c_1 \Omega^{-2} M_{1 1} \mu_{1 0} - c_1 \Omega^{-2} M_{1 1} - 2 c_1 \Omega^{-2} M_{2 2} - 2 c_1 \Omega^{-2} \mu_{1 0}^{2} M_{0 2} - c_1 \Omega^{-2} \mu_{0 1} \mu_{1 0} - 2 c_1 \Omega^{-2} M_{2 1} \mu_{0 1} - 4 c_1 \Omega^{-2} M_{1 2} \mu_{1 0} - 4 c_1 \Omega^{-2} M_{1 1} \mu_{0 1} \mu_{1 0}
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
latexify(closed_raw_eqs, :closure)
```
```math
\begin{align*}
\mu{_{21}} =& \mu{_{01}} \mu{_{20}} + 2 \mu{_{10}} \mu{_{11}} - 2 \mu{_{01}} \mu{_{10}}^{2} \\
\mu{_{12}} =& \mu{_{02}} \mu{_{10}} + 2 \mu{_{01}} \mu{_{11}} - 2 \mu{_{10}} \mu{_{01}}^{2} \\
\mu{_{31}} =& \mu{_{01}} \mu{_{30}} + 6 \mu{_{01}} \mu{_{10}}^{3} + 3 \mu{_{10}} \mu{_{21}} + 3 \mu{_{11}} \mu{_{20}} - 6 \mu{_{11}} \mu{_{10}}^{2} - 6 \mu{_{01}} \mu{_{10}} \mu{_{20}} \\
\mu{_{22}} =& \mu{_{02}} \mu{_{20}} + 2 \mu{_{11}}^{2} + 2 \mu{_{01}} \mu{_{21}} + 2 \mu{_{10}} \mu{_{12}} + 6 \mu{_{01}}^{2} \mu{_{10}}^{2} - 2 \mu{_{02}} \mu{_{10}}^{2} - 2 \mu{_{20}} \mu{_{01}}^{2} - 8 \mu{_{01}} \mu{_{10}} \mu{_{11}}
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
M{_{31}} =& 3 M{_{11}} M{_{20}} \\
M{_{22}} =& M{_{02}} M{_{20}} + 2 M{_{11}}^{2}
\end{align*}
```
Higher order central moments under normal closure take a rather simple form compared to their raw moment equivalents. That can be expected due to the relationship between central moments and cumulants, on which [the closure is based](@ref normal_closure).

## Solving Moment Equations

The closed moment equations can be solved numerically using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl/) that provides a variety of highly-efficient solvers and analysis tools. In order to do so, first we need to specify the values of all system parameters, the initial condition and the timespan to solve over. Then the [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) corresponding to the moment equations can be [directly converted](https://mtk.sciml.ai/stable/#Compatible-Numerical-Solvers-1) into an [`ODEProblem`](@ref) that can finally be solved. Let's go through the procedure step-by-step for the closed raw moment equations (`closed_raw_eqs`). Most of what is covered below is based on [on Catalyst's introductory tutorial](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/introduction_to_catalyst/#Stochastic-simulation-algorithms-(SSAs)-for-stochastic-chemical-kinetics).

We start with the parameters. Note that the internal parameter ordering [can vary](https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#Why-are-my-parameters-some-obscure-object?), and we need to build a mapping from the symbolic parameters to their respective numerical values. This can be done as:
```julia
pmap = [:c₁ => 0.9, :c₂ => 2, :c₃ => 1, :c₄ => 1, :Ω => 100]
```
The next ingredient, the time interval to solve on, can be specified simply as:
```julia
tspan = (0., 100.)
```
Lastly, we can specify the initial condition. Usually, when working with moment equations, we consider *deterministic* initial conditions so that the molecule numbers at initial time take the specified integer values with probability one. We can define the initial molecule numbers as $X(t=0) = X_0$ and $Y(t=0) = Y_0$. Probability one implies that initially the means will be equal to the molecule numbers, i.e., $μ_{10}(t=0) = X_0$ and $μ_{01}(t=0) = Y_0$, whereas all higher order raw moments will be products of the corresponding powers of the means, e.g., $μ_{21} = X_0^2 Y_0$. Such variable mapping under deterministic initial conditions is done automatically by default when the initial molecule numbers are passed to the [`ODEProblem`](@ref) constructor:
```julia
# initial molecule numbers
u0map = [:X => 1, :Y => 1]
# create the ODEProblem
oprob = ODEProblem(closed_raw_eqs, u0map, tspan, pmap)
```
Finally, we have everything we need to solve the raw moment equations which can be done using any ODE solver [implemented within DifferentialEquations.jl](https://docs.sciml.ai/Catalyst/stable/model_simulation/simulation_introduction/#simulation_intro_solver_options). We choose the commonly used `Tsit5()` solver and then [plot](https://diffeq.sciml.ai/stable/basics/plot/#plot) the obtained mean molecule numbers:
```julia
# loading only the solver-specific library (faster)
using OrdinaryDiffEqTsit5
sol = solve(oprob, Tsit5(), saveat=0.1)

using Plots
plot(sol, idxs=[1,2], lw=2)
```
![Brusselator means 1](../assets/brusselator_means_1.svg)

The obtained moment dynamics show damped oscillations which is the expected averaged behaviour of the Brusselator in a vast swathe of parameter space [1]. However, to establish more clearly how well a second order moment expansion with normal moment closure performs for this system and this specific set of parameters, we can compare the resulting moment trajectories to the true moment estimates obtained using [Gillespie's Stochastic Simulation Algorithm (SSA)](https://en.wikipedia.org/wiki/Gillespie_algorithm).

## Stochastic Simulation

To run the SSA for a given reaction network, we build a [JumpProcesses](https://github.com/SciML/JumpProcesses.jl)  [`JumpProblem`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#defining_jump_problem) using Gillespie's `Direct` method. Note that other SSA variants are also available, see the [documentation](https://diffeq.sciml.ai/dev/types/jump_types/#Constant-Rate-Jump-Aggregators). Moreover, in order to run many realisations of the jump process, we define a corresponding [`EnsembleProblem`](https://diffeq.sciml.ai/stable/features/ensemble/#ensemble). All of this can be done as follows:
```julia
using JumpProcesses

# convert ReactionSystem into JumpSystem
jsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
jsys = complete(jsys) 

# create a DiscreteProblem encoding that the molecule numbers are integer-valued
dprob = DiscreteProblem(jsys, u0map, tspan, pmap) # same parameters as defined earlier

# create a JumpProblem: specify Gillespie's Direct Method as the solver
# and SET save_positions to (false, false) as otherwise time of each
# reaction occurence would be saved (complicating moment estimates)
jprob = JumpProblem(jsys, dprob, Direct(), save_positions=(false, false))

# define an EnsembleProblem to simulate multiple trajectories
ensembleprob = EnsembleProblem(jprob)

# simulate 10000 SSA trajectories
@time sol_SSA = solve(ensembleprob, SSAStepper(), saveat=0.1, trajectories=10000)
```
Now we use the DifferentialEquations [ensemble statistics tools](https://diffeq.sciml.ai/stable/features/ensemble/#Analyzing-an-Ensemble-Experiment) to calculate the SSA mean values and plot them:
```julia
using DiffEqBase.EnsembleAnalysis

means_SSA = timeseries_steps_mean(sol_SSA)
plot!(means_SSA, lw=2, labels=["SSA μ₁₀" "SSA μ₀₁"], linestyle=:dash,
      linecolor=[1 2], background_color_legend=nothing, legend=:bottomright)
```
![Brusselator SSA 1](../assets/brusselator_ssa_1.svg)

The comparison to the SSA reveals that the second-order moment expansion using normal closure captures the correct qualitative behaviour of the Brusselator and provides reasonably accurate moment estimates given this particular parameter set. Note, however, that moment closure approximations can lead to unphysical results and suffer from numerical instabilities, please see the [Common Issues tutorial](common_issues.md) for more details.
