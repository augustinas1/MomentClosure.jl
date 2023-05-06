# [Parameter Estimation of Diffusion Processes](@id parameter_estimation_SDE)

In this tutorial, we look at the problem of model parameter identification of a diffusion process given timeseries data of its moments. Namely, we use moment closure approximations (MAs) to make the optimisation process more efficient by reducing the model evaluation time. All credit for this tutorial goes to Flemming Holtorf!

We consider a noisy variation of the [Lotka-Volterra model](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations), describing the interaction between a predator and prey species:
```math
\begin{align*}
    \begin{bmatrix} dx \\ dy \end{bmatrix} = 
    \begin{bmatrix} \gamma_1 x(t) - \gamma_2  x(t)  y(t)  \\ \gamma_4 x(t)  y(t) - \gamma_3 y(t) - \frac{1}{2} y(t) \end{bmatrix} \, dt 
    + \begin{bmatrix}  \gamma_5 x(t) \\ 0 \end{bmatrix} \, dW_t
\end{align*}
```
We can define this system of SDEs using ModelingToolkit as follows:
```julia
using ModelingToolkit

@variables t, x(t), y(t)
@parameters γ1, γ2, γ3, γ4, γ5
γ = [γ1, γ2, γ3, γ4, γ5]  
drift_eqs = [Differential(t)(x) ~ γ[1] * x - γ[2] * x * y ;
             Differential(t)(y) ~ γ[4] * x * y - γ[3] * y - y*0.5]
diff_eqs = [γ[5]*x; 0]
LV = SDESystem(drift_eqs, diff_eqs, t, [x,y], γ, name = :LV)
```
Next, we generate some data used for the parameter estimation. Namely, we collect timeseries data of means and variances of both species:
```julia
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis

N_samples = 1000
Tf = 10
t_data = 0:0.2:Tf
p_true = [γ[1] => 1, γ[2] => 2, γ[3] => 1, γ[4] => 2, γ[5] => 0.1]
u0 = [1.0, 0.25]
LV_data = solve(EnsembleProblem(SDEProblem(LV, u0, (0.0, Tf), p_true)), saveat = t_data, trajectories = N_samples)
means, vars = timeseries_steps_meanvar(LV_data)
```
Now we are ready to test if we can estimate the model parameters $\gamma_1, \dotsc, \gamma_5$ solely from the data collected above. We first approach this parameter identification problem using the asymptotically exact approach of estimating the means and variances of the process with ensemble averages. Accordingly, we construct the following loss function:
```julia
using LinearAlgebra

LV_sde = SDEProblem(LV, u0, (0.0, Tf), zeros(5))
function obj(p)
    prob = remake(LV_sde, p = p)
    sol = solve(EnsembleProblem(prob), saveat = t_data, trajectories = 1000)
    sol_mean, sol_vars = timeseries_steps_meanvar(sol)
    obj = sum(norm(sol_mean[i] - means[i])^2 for i in 1:length(t_data))
    obj += 1e4*sum(norm(sol_vars[i] - vars[i])^2 for i in 1:length(t_data))
    return obj
end
```
We can use this loss function with any suitable optimisation routine to identify a reasonable choice of the model parameters. For example, we can use a very simple derivative-free optimizer (Nelder-Mead method) implemented in the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package. Since a single evaluation of the objective function requires sampling and hence is relatively expensive, we impose the constraint that the optimiser cannot run more than 2 minutes:
```julia
using Optim

pmap = Dict(i => findfirst(isequal(γ[i]), parameters(LV)) for i in 1:5)
p_init = Dict(γ[1] => 1.3, γ[2] => 1.5, γ[3] => 1.4, γ[4] => 2.2, γ[5] => 0.1)
p = [p_init[γ[pmap[i]]] for i in 1:5]
opt_sampling = Optim.optimize(obj, p, Optim.Options(time_limit = 120))
```
```julia
* Status: failure (exceeded time limit of 120.0)

 * Candidate solution
    Final objective value:     3.550929e+00

 * Found with
    Algorithm:     Nelder-Mead

 * Convergence measures
    √(Σ(yᵢ-ȳ)²)/n ≰ 1.0e-08

 * Work counters
    Seconds run:   121  (vs limit 120)
    Iterations:    139
    f(x) calls:    337
```
We can visualise the moment statistics obtained using the estimated parameters and compare to the true data:
```julia
using Plots

t_detail = collect(0:0.01:Tf) 
opt_sol = solve(EnsembleProblem(SDEProblem(LV, u0, (0.0, Tf), opt_sampling.minimizer)), saveat = t_detail, trajectories = 1000)
opt_means, opt_vars = timeseries_steps_meanvar(opt_sol)

mean_comp = scatter(t_data, [m[1] for m in means], color = :blue,
                    xlabel = "time", ylabel = "population size", 
                    grid = false, title = "means", label = "⟨x⟩ data")
scatter!(mean_comp, t_data, [m[2] for m in means], color = :red, label = "⟨y⟩ data")
plot!(mean_comp, t_detail, [m[1] for m in opt_means], linewidth = 2, color = :blue, label = "⟨x⟩ SDE model")
plot!(mean_comp, t_detail, [m[2] for m in opt_means], linewidth = 2, color = :red, label = "⟨y⟩ SDE model")

var_comp = scatter(t_data, [v[1] for v in vars], color = :blue, grid = false,
                   xlabel = "time", title = "variances", label = "σ²(x) data", legend = :topleft)
scatter!(var_comp, t_data, [v[2] for v in vars], color = :red, label = "σ²(y) data")
plot!(var_comp, t_detail, [v[1] for v in opt_vars], color = :blue, label = "σ²(x) SDE model")
plot!(var_comp, t_detail, [v[2] for v in opt_vars], color = :red, label = "σ²(y) SDE model")

plot(mean_comp, var_comp, size = (1200.0, 400.0))
```
![LV SDE fig1](../assets/LV_SDE_fig1.svg)

We observe that the identified parameters match the data reasonably well. However, one may suspect that the optimiser is converging to a local minimum as the fit is not perfect.

Now we approach the same model identification problem via MAs in the hope of cutting down model evaluation cost, allowing us to identify better parameters in the same (or less) time. To that end, we construct an approximation of the moment dynamics of the process assuming that the distribution of the system state is approximately log-normal over the simulation horizon (using [log-normal MA](@ref log-normal_closure)). Then we can implement a simple loss function by comparing the moments predicted by the approximate model with those obtained from data:
```julia
using MomentClosure

LV_moments = moment_closure(generate_raw_moment_eqs(LV, 2), "log-normal")
u0map = deterministic_IC(u0, LV_moments)
closed_moment_prob = ODEProblem(LV_moments, u0map, (0.0, Tf), zeros(5))
function obj_MCA(p)
    prob = remake(closed_moment_prob; p=p)
    sol = solve(prob, Tsit5(), saveat = t_data)
    if sol.retcode == ReturnCode.Success
        obj = sum(norm(sol.u[i][1:2] - means[i])^2 for i in 1:length(t_data))
        obj += 1e4*sum((sol.u[i][3] - sol.u[i][1]^2  - vars[i][1])^2 for i in 1:length(t_data))
        obj += 1e4*sum((sol.u[i][5] - sol.u[i][2]^2  - vars[i][2])^2 for i in 1:length(t_data))
    else
        obj = 1e6
    end
    return obj
end
```
As before, any suitable optimisation routine can now be used to identify parameter values that result in a match between data and model prediction:
```julia
pmap = Dict(i => parameters(LV_moments.odes)[i] for i in 1:5)
p_init = Dict(γ[1] => 1.3, γ[2] => 1.5, γ[3] => 1.4, γ[4] => 2.2, γ[5] => 0.1)
p = [p_init[pmap[i]] for i in 1:5]
opt_MCA = Optim.optimize(obj_MCA, p, Optim.Options(time_limit = min(120, opt_sampling.time_run)))
```
```julia
* Status: success

 * Candidate solution
    Final objective value:     1.160452e+00

 * Found with
    Algorithm:     Nelder-Mead

 * Convergence measures
    √(Σ(yᵢ-ȳ)²)/n ≤ 1.0e-08

 * Work counters
    Seconds run:   0  (vs limit 120)
    Iterations:    362
    f(x) calls:    751
```
By visualising the results we see that indeed the identified parameters now provide a better match between data and the model, even if the original SDE model is being evaluated (note that the moment equations using log-normal MA provide a reasonably accurate approximation to the ensemble averages):
```julia
p_opt = [pmap[i] => opt_MCA.minimizer[i] for i in 1:5]
t_detail = collect(0:0.01:Tf) 
opt_sol = solve(EnsembleProblem(SDEProblem(LV, u0, (0.0, Tf), p_opt)), saveat = t_detail, trajectories = 1000)
opt_means = [timestep_mean(opt_sol, i) for i in 1:length(t_detail)]
opt_vars = [timestep_meanvar(opt_sol, i)[2] for i in 1:length(t_detail)]

opt_sol_approx = solve(ODEProblem(LV_moments, u0map, (0.0, Tf), p_opt), saveat = t_detail)

mean_comp = scatter(t_data, [m[1] for m in means], color = :blue,
                     xlabel = "time", ylabel = "population size", 
                     grid = false, title = "means", label = "⟨x⟩ data")
scatter!(mean_comp, t_data, [m[2] for m in means], color = :red, label = "⟨y⟩ data")
plot!(mean_comp, t_detail, [m[1] for m in opt_means], linewidth = 2, color = :blue, label = "⟨x⟩ SDE model")
plot!(mean_comp, t_detail, [m[2] for m in opt_means], linewidth = 2, color = :red, label = "⟨y⟩ SDE model")
plot!(mean_comp, t_detail, [m[1] for m in opt_sol_approx.u], linewidth = 2, color = :black, linestyle = :dash, label = "log-normal MA")
plot!(mean_comp, t_detail, [m[2] for m in opt_sol_approx.u], linewidth = 2, color = :black, linestyle = :dash, label = nothing)

var_comp = scatter(t_data, [v[1] for v in vars], color = :blue,
                   xlabel = "time", title = "variances", grid = false, label = "σ²(x) data", legend = :topleft)
scatter!(var_comp, t_data, [v[2] for v in vars], color = :red, label = "σ²(y) data")
plot!(var_comp, t_detail, [v[1] for v in opt_vars], color = :blue, label = "σ²(x) SDE model")
plot!(var_comp, t_detail, [v[2] for v in opt_vars], color = :red, label = "σ²(y) SDE model")
plot!(var_comp, t_detail, [m[3] - m[1]^2 for m in opt_sol_approx.u], linewidth = 2, color = :black, linestyle = :dash, label = "log-normal MA")
plot!(var_comp, t_detail, [m[5] - m[2]^2 for m in opt_sol_approx.u], linewidth = 2, color = :black, linestyle = :dash, label = nothing)

plot(mean_comp, var_comp, size = (1200.0, 400.0))
```
![LV SDE fig2](../assets/LV_SDE_fig2.svg)

Finally, note that the parameter estimation is much faster using moment equations (compare to 120 s using the basic approach)
```julia
opt_MCA.time_run
```
```julia
0.10686397552490234
```