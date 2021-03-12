# [Common Issues](@id common_issues)

It is important to be well aware that moment closure approximations are based on *ad hoc* assumptions and no rigorous and general predictions can be made on whether the results will be accurate or even physically meaningful [1-2]. Moreover, the truncated moment equations are prone to numerical instabilities and it may not be possible to solve them for the entire time course [3]. In this tutorial, we walk through a number of such issues encountered in the analysis of the Brusselator introduced in the [previous tutorial](using_momentclosure.md).

We first redefine the system and its parameters for completeness:
```julia
using MomentClosure, Catalyst, OrdinaryDiffEq, Plots

rn = @reaction_network begin
  (c₁/Ω^2), 2X + Y → 3X
  (c₂), X → Y
  (c₃*Ω, c₄), 0 ↔ X
end c₁ c₂ c₃ c₄ Ω

p = [0.9, 2, 1, 1, 100]
u₀ = [1, 1]
tspan = (0., 100.)

raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=false)
```
As we have seen earlier, second-order moment expansion using normal closure approximates the true system dynamics sufficiently accurately but it would be interesting to see how other closures compare. Let's try applying zero closure:
```julia
closed_raw_eqs = moment_closure(raw_eqs, "zero")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Tsit5(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2)
```
![Brusselator issue 1](../assets/brusselator_issue_1.svg)

The trajectory of $μ₀₁$ becomes negative and so zero closure fails to provide physically meaningful results for this parameter set. Note that is important to correctly choose the ODE solver depending on the stiffness of the system and the accuracy required. We tried a number of [recommended DifferentialEquations solvers](https://diffeq.sciml.ai/stable/solvers/ode_solve/) here but none seemed to improve the results.

We apply log-normal closure next and try solving the equations using `Tsit5()` solver which immediately throws a `DomainError`: the mean trajectories become negative again and they are passed to $\log$ terms in the higher order moment closure functions which result in the error message. Nevertheless, this can be overcome by using one of the solvers advised for stiff problems, e.g., `Rodas4P()`:
```julia
closed_raw_eqs = moment_closure(raw_eqs, "log-normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Rodas4P(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2, legend=:bottomright)
```
![Brusselator issue 2](../assets/brusselator_issue_2.svg)

We observe sustained oscillatory behaviour instead of the expected damped oscillations. This result is unphysical: single SSA trajectories (that may display sustained oscillations) get dephased over time and hence the ensemble average should always show damped or overdamped oscillations [1].

Normal closure is also quite fragile. This can be seen by simply including the combinatorial scaling of the mass-action propensity functions with `combinatoric_ratelaw=true` which leads to unphysical sustained oscillatory trajectories:
```julia
raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=true)
closed_raw_eqs = moment_closure(raw_eqs, "normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Tsit5(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2)
```
![Brusselator issue 3](../assets/brusselator_issue_3.svg)

Nevertheless, this can be improved upon by increasing the order of moment expansion:
```julia
raw_eqs = generate_raw_moment_eqs(rn, 3, combinatoric_ratelaw=true)
closed_raw_eqs = moment_closure(raw_eqs, "normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Tsit5(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2, legend=:bottomright)
```
![Brusselator issue 4](../assets/brusselator_issue_4.svg)

Some dampening in the system is now visible. Increasing the expansion order to `4` finally leads to sensible results:
```julia
raw_eqs = generate_raw_moment_eqs(rn, 4, combinatoric_ratelaw=true)
closed_raw_eqs = moment_closure(raw_eqs, "normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Tsit5(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2)
```
![Brusselator issue 5](../assets/brusselator_issue_5.svg)

For dessert, we consider unphysical divergent trajectories—a frequent problem with moment equations [3]. A good example is the second-order moment expansion including the combinatorial scaling of propensities with log-normal closure applied:
```julia
raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=true)
closed_raw_eqs = moment_closure(raw_eqs, "log-normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Rodas4P(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2)
```
![Brusselator issue 6](../assets/brusselator_issue_6.svg)

In contrast to normal closure, increasing the expansion order makes the problem worse:
```julia
raw_eqs = generate_raw_moment_eqs(rn, 3, combinatoric_ratelaw=true)
closed_raw_eqs = moment_closure(raw_eqs, "log-normal")

u₀map = deterministic_IC(u₀, closed_raw_eqs)
oprob = ODEProblem(closed_raw_eqs, u₀map, tspan, p)
sol = solve(oprob, Rodas4P(), saveat=0.1)

plot(sol, vars=(0, [1,2]), lw=2)
```
```julia
┌ Warning: Interrupted. Larger maxiters is needed.
└ @ SciMLBase C:\Users\asukys\.julia\packages\SciMLBase\Afx1r\src\integrator_interface.jl:331
```
![Brusselator issue 7](../assets/brusselator_issue_7.svg)

Note that the solver throws a warning being unable to evaluate the trajectories for the entire time course (other solvers perform similarly in this case). This usually implies that the moment ODE system is too stiff and cannot be solved: the time derivatives grow unboundedly and the solver timestep is being constantly reduced, requiring an ever-increasing number of the solver iterations (hence the [`maxiters`](https://diffeq.sciml.ai/stable/basics/common_solver_opts/#Miscellaneous) warning).

## References

[1]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Comparison of different moment-closure approximations for stochastic chemical kinetics", The Journal of Chemical Physics 143, 185101 (2015). [https://doi.org/10.1063/1.4934990](https://doi.org/10.1063/1.4934990)

[2]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Validity conditions for moment closure approximations in stochastic chemical kinetics", The Journal of Chemical Physics 141, 084103 (2014). [https://doi.org/10.1063/1.4892838](https://doi.org/10.1063/1.4892838)

[3]: E. Lakatos, A. Ale, P. D. W. Kirk, and M. P. H. Stumpf, "Multivariate moment closure techniques for stochastic kinetic models", The Journal of Chemical Physics 143, 094107 (2015). [https://doi.org/10.1063/1.4929837](https://doi.org/10.1063/1.4929837)
