# [Geometrically Distributed Reaction Products and Conditional Closures](@id geometric-and-conditional)

In this tutorial, we demonstrate using [Catalyst](https://github.com/SciML/Catalyst.jl) how to define chemical systems involving reactions which products are geometrically distributed random variables. As an example, we consider an autoregulatory (repressive) genetic feedback loop where proteins are expressed in bursts with a geometric burst size distribution, as described by Soltani et al. [1]. Moreover, the state of the gene is modelled as a binary variable—we demonstrate how conditional derivative matching and conditional gaussian closures can be used to approximate such systems, in turn reproducing a number of results from Ref. [1].

We consider a negative feedback loop described by the following reactions:
```math
\begin{align*}
G^* &\stackrel{k_{on}}{\rightarrow} G, \\
G + &\stackrel{k_{off}*P^2}{\rightarrow} G^*, \\
G  &\stackrel{k_p}{\rightarrow} G + mP, \\
P &\stackrel{\gamma_p}{\rightarrow} ∅.
\end{align*}
```
A gene in the network switches between ON ($G$) and OFF ($G^{\\*}$) states: proteins are produced in the transcriptionally active ON state but the gene can be turned OFF by two protein molecules binding to the promoter region and thus blocking transcription (proteins decay at a constant rate irrespective of the gene state). Note that the gene state can be interpreted as a distinct species that have either zero or one copy number per cell. In other words, it is a Bernoulli random variable: $0$ in the OFF state and $1$ in the ON state.

The transcription (mRNA) dynamics are not modelled explicitly in this gene circuit. Instead, under the assumption of fast mRNA decay, proteins are taken to be produced in bursts of size $m$, where $m$ is a random variable sampled from the geometric distribution $\phi(m) = b^m/(1+b)^{m+1}$ (with $b$ denoting the mean burst size) [2]. Chemical reaction networks that include reactions which products are independent geometrically distributed random variables (or other symbolic variables) can be defined using [Catalyst](https://github.com/SciML/Catalyst.jl) as demonstrated [in this tutorial](https://docs.sciml.ai/Catalyst/stable/model_creation/parametric_stoichiometry/#Gene-expression-with-randomly-produced-amounts-of-protein). Note that the referenced tutorial refers back to this tutorial because in the ancient times Catalyst did not provided such functionality and the previous iteration of MomentClosure implemented a `ReactionSystemMod` type that offered limited support for such systems (now deprecated as Catalyst does it better in a unified API).

Using `Catalyst`, our gene network model can be constructed as follows:
```julia
# load all the packages we will need
using MomentClosure, Catalyst, Distributions, JumpProcesses, DiffEqBase, OrdinaryDiffEq, DiffEqBase.EnsembleAnalysis, Plots, Latexify

# Proteins are produced in bursts of size m, where m is a geometric random variable with mean b.
# Note that if b is the mean burst size, then p = 1/(1+b). 
# Implemented by first registering the distribution with Symbolics
@register_symbolic Distributions.Geometric(b)
@parameters b
m = rand(Distributions.Geometric(1/(1+b)))

rn = @reaction_network begin
      @parameters k_on k_off k_p γ_p
      k_on*(1-g), 0 --> g  # G* -> G
      k_off*P^2, g --> 0   # G -> G*
      k_p, g --> g + $m*P  # G -> G + mP, m ~ Geometric(p)
      γ_p, P --> 0         # P -> ∅
end
```
We can now generate the raw moment equations up to third order:
```julia
eqs = generate_raw_moment_eqs(rn, 3)
latexify(eqs)
```
```math
\begin{align*}
\frac{d\mu_{1 0}}{dt} =& k_{on} - k_{off} \mu_{1 2} - k_{on} \mu_{1 0} \\
\frac{d\mu_{0 1}}{dt} =& b k_{p} \mu_{1 0} - \gamma_{p} \mu_{0 1} \\
\frac{d\mu_{2 0}}{dt} =& k_{on} + k_{off} \mu_{1 2} + k_{on} \mu_{1 0} - 2 k_{on} \mu_{2 0} - 2 k_{off} \mu_{2 2} \\
\frac{d\mu_{1 1}}{dt} =& k_{on} \mu_{0 1} + b k_{p} \mu_{2 0} - k_{on} \mu_{1 1} - \gamma_{p} \mu_{1 1} - k_{off} \mu_{1 3} \\
\frac{d\mu_{0 2}}{dt} =& \gamma_{p} \mu_{0 1} + b k_{p} \mu_{1 0} + 2 b k_{p} \mu_{1 1} + 2 k_{p} b^{2} \mu_{1 0} - 2 \gamma_{p} \mu_{0 2} \\
\frac{d\mu_{3 0}}{dt} =& k_{on} + 3 k_{off} \mu_{2 2} + 2 k_{on} \mu_{1 0} - k_{off} \mu_{1 2} - 3 k_{off} \mu_{3 2} - 3 k_{on} \mu_{3 0} \\
\frac{d\mu_{2 1}}{dt} =& k_{on} \mu_{0 1} + k_{on} \mu_{1 1} + k_{off} \mu_{1 3} + b k_{p} \mu_{3 0} - 2 k_{on} \mu_{2 1} - \gamma_{p} \mu_{2 1} - 2 k_{off} \mu_{2 3} \\
\frac{d\mu_{1 2}}{dt} =& k_{on} \mu_{0 2} + \gamma_{p} \mu_{1 1} + b k_{p} \mu_{2 0} + 2 k_{p} b^{2} \mu_{2 0} + 2 b k_{p} \mu_{2 1} - k_{off} \mu_{1 4} - k_{on} \mu_{1 2} - 2 \gamma_{p} \mu_{1 2} \\
\frac{d\mu_{0 3}}{dt} =& b k_{p} \mu_{1 0} + 3 \gamma_{p} \mu_{0 2} + 3 b k_{p} \mu_{1 1} + 3 b k_{p} \mu_{1 2} + 6 k_{p} b^{2} \mu_{1 0} + 6 k_{p} b^{3} \mu_{1 0} + 6 k_{p} b^{2} \mu_{1 1} - \gamma_{p} \mu_{0 1} - 3 \gamma_{p} \mu_{0 3}
\end{align*}
```
A lot of information in this system of ODEs is redundant as the gene state is a Bernoulli variable that (in our case) has the following properties:
```math
\begin{align*}
μ_{j0} &= μ_{10}, \quad j \in {2, 3, \dotsc}, \\
μ_{jk} &= μ_{1k}, \quad j,k \in {1, 2, 3, \dotsc}
\end{align*}
```
We can *clean* the moment equations accordingly by specifying which molecule numbers are Bernoulli variables and using the [`bernoulli_moment_eqs`](@ref) function:
```julia
# specify the indices of species which molecule numbers are binary
binary_vars = [1]
# simplify the moment equations using properties of Bernoulli variables
clean_eqs = bernoulli_moment_eqs(eqs, binary_vars)
latexify(clean_eqs)
```
```math
\begin{align*}
\frac{d\mu_{1 0}}{dt} =& k_{on} - k_{off} \mu_{1 2} - k_{on} \mu_{1 0} \\
\frac{d\mu_{0 1}}{dt} =& b k_{p} \mu_{1 0} - \gamma_{p} \mu_{0 1} \\
\frac{d\mu_{1 1}}{dt} =& k_{on} \mu_{0 1} + b k_{p} \mu_{1 0} - k_{on} \mu_{1 1} - \gamma_{p} \mu_{1 1} - k_{off} \mu_{1 3} \\
\frac{d\mu_{0 2}}{dt} =& \gamma_{p} \mu_{0 1} + b k_{p} \mu_{1 0} + 2 b k_{p} \mu_{1 1} + 2 k_{p} b^{2} \mu_{1 0} - 2 \gamma_{p} \mu_{0 2} \\
\frac{d\mu_{1 2}}{dt} =& k_{on} \mu_{0 2} + \gamma_{p} \mu_{1 1} + b k_{p} \mu_{1 0} + 2 b k_{p} \mu_{1 1} + 2 k_{p} b^{2} \mu_{1 0} - k_{off} \mu_{1 4} - k_{on} \mu_{1 2} - 2 \gamma_{p} \mu_{1 2} \\
\frac{d\mu_{0 3}}{dt} =& b k_{p} \mu_{1 0} + 3 \gamma_{p} \mu_{0 2} + 3 b k_{p} \mu_{1 1} + 3 b k_{p} \mu_{1 2} + 6 k_{p} b^{2} \mu_{1 0} + 6 k_{p} b^{3} \mu_{1 0} + 6 k_{p} b^{2} \mu_{1 1} - \gamma_{p} \mu_{0 1} - 3 \gamma_{p} \mu_{0 3}
\end{align*}
```
The system of ODEs is now much simpler and we can see that there are two higher-order moments we need to truncate: $\mu_{13}$ and $μ_{14}$. We consider normal, derivative matching, conditional gaussian and conditional derivative matching closures to see how well they compare. First we apply different closures and print out the corresponding higher-order moment expressions in order to check that our results are consistent with those published in [1].

For normal closure:
```julia
normal_eqs = moment_closure(eqs, "normal", binary_vars)
latexify(normal_eqs, :closure)
```
```math
\begin{align*}
\mu{_{13}} =& \mu{_{03}} \mu{_{10}} + 3 \mu{_{01}} \mu{_{12}} + 3 \mu{_{02}} \mu{_{11}} + 6 \mu{_{10}} \mu{_{01}}^{3} - 6 \mu{_{11}} \mu{_{01}}^{2} - 6 \mu{_{01}} \mu{_{02}} \mu{_{10}} \\
\mu{_{14}} =& \mu{_{04}} \mu{_{10}} + 4 \mu{_{01}} \mu{_{13}} + 6 \mu{_{02}} \mu{_{12}} + 4 \mu{_{03}} \mu{_{11}} + 24 \mu{_{11}} \mu{_{01}}^{3} + 36 \mu{_{02}} \mu{_{10}} \mu{_{01}}^{2} - 24 \mu{_{10}} \mu{_{01}}^{4} - 6 \mu{_{10}} \mu{_{02}}^{2} - 12 \mu{_{12}} \mu{_{01}}^{2} - 24 \mu{_{01}} \mu{_{02}} \mu{_{11}} - 8 \mu{_{01}} \mu{_{03}} \mu{_{10}}
\end{align*}
```
Note that our expression for $μ_{14}$ is different from the equivalent $\langle gp^4 \rangle$ in Eq. (8) of [1] as it seems that the authors introduced a mistake considering the central moment and cumulant relationship in Eq. (7).

Using derivative matching, we obtain:
```julia
dm_eqs = moment_closure(eqs, "derivative matching", binary_vars)
latexify(dm_eqs, :closure)
```
```math
\begin{align*}
\mu_{1 3} =& \frac{\mu_{0 1}^{3} \mu_{1 2}^{3} \mu_{0 3} \mu_{1 0}}{\mu_{0 2}^{3} \mu_{1 1}^{3}} \\
\mu_{1 4} =& \frac{\mu_{0 2}^{6} \mu_{1 1}^{4} \mu_{1 3}^{4} \mu_{0 4}}{\mu_{0 1}^{4} \mu_{0 3}^{4} \mu_{1 2}^{6} \mu_{1 0}}
\end{align*}
```
For conditional gaussian closure:
```julia
cond_gaussian_eqs = moment_closure(eqs, "conditional gaussian", binary_vars)
latexify(cond_gaussian_eqs, :closure)
```
```math
\begin{align*}
\mu_{1 3} =& \frac{-2 \mu_{1 1}^{3}}{\mu_{1 0}^{2}} + \frac{3 \mu_{1 1} \mu_{1 2}}{\mu_{1 0}} \\
\mu_{1 4} =& \frac{3 \mu_{1 2}^{2}}{\mu_{1 0}} + \frac{6 \mu_{1 1}^{4}}{\mu_{1 0}^{3}} + \frac{-12 \mu_{1 1}^{2} \mu_{1 2}}{\mu_{1 0}^{2}} + \frac{4 \mu_{1 1} \mu_{1 3}}{\mu_{1 0}}
\end{align*}
```
And, finally, for conditional derivative matching:
```julia
cond_dm_eqs = moment_closure(eqs, "conditional derivative matching", binary_vars)
latexify(cond_dm_eqs, :closure)
```
```math
\begin{align*}
\mu_{1 3} =& \frac{\mu_{0 1}^{3} \mu_{1 2}^{3} \mu_{0 3} \mu_{1 0}}{\mu_{0 2}^{3} \mu_{1 1}^{3}} \\
\mu_{1 4} =& \frac{\mu_{0 2}^{6} \mu_{1 1}^{4} \mu_{1 3}^{4} \mu_{0 4}}{\mu_{0 1}^{4} \mu_{0 3}^{4} \mu_{1 2}^{6} \mu_{1 0}}
\end{align*}
```
All these results are consistent with [1], reassuring that the model and closures are implemented as intended. Finally, we can proceed to solve the resulting ODEs and compare the resulting means and standard deviations to the SSA. Following Soltani et al. [1], we define all model parameters as:
```julia
mean_p = 200
mean_b = 70
γ_p_val = 1
k_off_val = 0.001
k_on_val = 0.05

k_p_val = mean_p * γ_p_val * (k_off_val * mean_p^2 + k_on_val) / (k_on_val * mean_b)

pmap = [:k_on => k_on_val,
          :k_off => k_off_val,
          :k_p => k_p_val,
          :γ_p => γ_p_val,
          :b => mean_b]

# initial gene state and protein number, order [g, p]
u₀ = [1, 1]

# time interval to solve on
tspan = (0., 6.0)
```
The reaction network with geometric bursts can be simulated using the SSA as usual:
```julia
# convert the reaction network into a system of jump processes
jsys = convert(JumpSystem, rn; combinatoric_ratelaws=false)
jsys = complete(jsys)

# create a discrete problem setting the simulation parameters
dprob = DiscreteProblem(u₀, tspan, pmap)

# create a JumpProblem compatible with ReactionSystemMod
jprob = JumpProblem(rn, dprob, Direct(), save_positions=(false, false))

# simulate 2×10⁴ SSA trajectories
ensembleprob  = EnsembleProblem(jprob)
@time sol_SSA = solve(ensembleprob, SSAStepper(), saveat=0.1, trajectories=20000)
# compute the means and variances
means_ssa, vars_ssa = timeseries_steps_meanvar(sol_SSA)
```
```julia
4.858357 seconds (14.30 M allocations: 901.032 MiB, 6.25% gc time, 63.45% compilation time)
```
We continue to solve the moment equations for each closure:
```julia
plt_m = plot()   # plot mean protein number
plt_std = plot() # plot ssd of protein number

# construct the initial molecule number mapping
u₀map = deterministic_IC(u₀, dm_eqs)

# solve moment ODEs for each closure and plot the results
for closure in ["normal", "derivative matching",
                "conditional gaussian", "conditional derivative matching"]

    # it is very quick so we just apply all closures again
    closed_eqs = moment_closure(eqs, closure, binary_vars)

    # solve the system of moment ODEs
    oprob = ODEProblem(closed_eqs, u₀map, tspan, pmap)
    sol = solve(oprob, AutoTsit5(Rosenbrock23()), saveat=0.01)

    # μ₀₁ is 2nd and μ₀₂ is 4th element in sol
    # can check the order with `closed_eqs.odes.states`
    plt_m = plot!(plt_m, sol, idxs=[2], label=closure)
    plt_std = plot!(plt_std, sol.t, sqrt.(sol[4, :] .- sol[2, :].^2), label=closure)

end

plt_m = plot!(plt_m, xlabel="Time [hr]", ylabel="Protein mean level")
plt_m = plot!(plt_m, means_ssa.t, means_ssa[2,:], label="SSA", linestyle=:dash, color="gray")
plt_std = plot!(plt_std, xlabel="Time [hr]", ylabel="Protein standard deviation")
plt_std = plot!(plt_std, vars_ssa.t, sqrt.(vars_ssa[2,:]), label="SSA", linestyle=:dash, color="gray")
```
```julia
plot(plt_m, lw=2)
```
![Gene 1 means](../assets/gene_1_means.svg)

```julia
plot(plt_std, lw=2)
```
![Gene 2 stds](../assets/gene_1_stds.svg)

We observe that the two conditional closures give the most accurate results. Therefore, closing moments conditioned on the gene being active is an effective approximation in this case. Note that our results for normal closure are different from [1], most likely due to the incorrect fifth-order moment closure function given in the paper.

We now turn to a more complex network involving two genes (two Bernoulli variables), a repressor-activator circuit: an activator protein $Y$ expressed by one gene can turn ON another gene which produces a repressor protein $X$, that subsequently can bind to the gene promoter region of the activator-producing gene turning it OFF (see [1] for more details). As before, $X$ and $Y$ proteins are expressed in geometrically-distributed bursts. The system can be specified as follows:
```julia
@parameters b_x b_y
# could redefine in terms of success probabilities p_x and p_y as above
m = rand(Distributions.Geometric(1/(1+b_x)))
l = rand(Distributions.Geometric(1/(1+b_y)))

# g_x - gene state of X protein producing gene
# g_y - gene state of Y protein producing gene
# x, y - proteins X and Y
rn = @reaction_network begin
      @parameters kx_on kx_off ky_on ky_off k_x γ_x k_y γ_y
      kx_on*(1-g_x)*y, 0 --> g_x  # 0   -> g_x
      kx_off,          g_x --> 0  # g_x -> 0
      ky_on*(1-g_y),   0 --> g_y  # 0 -> g_y
      ky_off*x,        g_y --> 0  # g_y -> 0
      k_x*g_x,         0 --> $m*x # 0 -> mx, m ~ Geometric(mean=b_x)
      γ_x,             x --> 0    # x -> 0
      k_y*g_y,         0 --> $l*y # 0 -> ly, l ~ Geometric(mean_b_y)
      γ_y,             y --> 0    # y -> 0
end

# both g_x and g_y are Bernoulli random variables
binary_vars = [1, 2]

# Parameter initialisation

mean_x = 100
mean_y = 100
mean_b_x = 5
mean_b_y = 5
γ_x_val = 1
γ_y_val = 1
kx_off_val = 4
ky_on_val = 0.3
kx_on_val = 0.05
ky_off_val = 0.05

k_x_val = mean_x * γ_x_val * (kx_off_val * mean_y^2 + kx_on_val) / (kx_on_val * mean_b_x)
k_y_val = mean_y * γ_y_val * (ky_off_val * mean_x^2 + ky_on_val) / (ky_on_val * mean_b_y)

# unclear if Soltani et al. (2015) actually used this parameter set as X numbers
# jump to millions making SSA extremely slow...

# introduce additional rescaling (otherwise rate coefficients are too high)
k_x_val *= 0.00003
k_y_val *= 0.01

# parameter mapping
pmap = [:kx_on => kx_on_val,
          :kx_off => kx_off_val,
          :ky_on => ky_on_val,
          :ky_off => ky_off_val,
          :k_x => k_x_val,
          :k_y => k_y_val,
          :γ_x => γ_x_val,
          :γ_y => γ_y_val,
          :b_x => mean_b_x,
          :b_y => mean_b_y]

# initial gene state and protein number, order [g_x, g_y, x, y]
u₀ = [1, 1, 1, 1]

# time interval to solve on
tspan = (0., 12.0)
```
Note that here we use a different parameter set from the one considered in [1] as otherwise the number of $X$ proteins goes into millions, making the comparison between closure results difficult due to SSA being extremely slow—it appears that either we have included a mistake setting the parameter values or there is a typo in the paper.

We can run SSA as follows:
```julia
jsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
jsys = complete(jsys)
dprob = DiscreteProblem(jsys, u₀, tspan, pmap)
jprob = JumpProblem(jsys, dprob, Direct(), save_positions=(false, false))

ensembleprob  = EnsembleProblem(jprob)
@time sol_SSA = solve(ensembleprob, SSAStepper(), saveat=0.1, trajectories=10000)
means_ssa, vars_ssa = timeseries_steps_meanvar(sol_SSA)
```
```julia
13.712331 seconds (43.27 M allocations: 1.697 GiB, 6.56% gc time, 10.01% compilation time)
```
We use fourth order moment expansion and apply different closure methods as done for the negative feedback loop, considering the number of activator proteins $Y$ and its standard deviation:
```julia
plt_m = plot()   # plot mean activator protein number
plt_std = plot() # plot ssd of activator protein number

for closure in ["derivative matching", "conditional derivative matching"]

    closed_eqs = moment_closure(eqs, closure, binary_vars)

    u₀map = deterministic_IC(u₀, closed_eqs)
    oprob = ODEProblem(closed_eqs, u₀map, tspan, pmap)
    sol = solve(oprob, Tsit5(), saveat=0.1)

    # μ₀₀₀₁ is the 4th and μ₀₀₀₂ is the 12th element in sol (can check with closed_eqs.odes.states)
    plt_m = plot!(plt_m, sol, idxs=[4], label=closure)
    plt_std = plot!(plt_std, sol.t, sqrt.(sol[12, :] .- sol[4, :].^2), label=closure)
end

plt_m = plot!(plt_m, xlabel="Time [hr]", ylabel="Activator mean level")
plt_m = plot!(plt_m, means_ssa.t, means_ssa[4,:], label="SSA", linestyle=:dash, color="gray")
plt_std = plot!(plt_std, xlabel="Time [hr]", ylabel="Activator standard deviation")
plt_std = plot!(plt_std, vars_ssa.t, sqrt.(vars_ssa[4,:]), label="SSA", linestyle=:dash, color="gray")
```
Note that here we apply only derivative matching and conditional derivative matching, as both normal and conditional gaussian closures failed to be evaluated for the entire time course (irrespective of the solver used). Finally, we can plot the results:

```julia
plot(plt_m, lw=2)
```
![Gene 2 means](../assets/gene_2_means.svg)

```julia
plot(plt_std, lw=2, xlims=(0., 12.))
```
![Gene 2 stds](../assets/gene_2_stds.svg)

The trajectories obtained with conditional derivative matching quite accurately match the true dynamics, showing that it is indeed a powerful closure method for reaction networks involving species with very low (binary) copy numbers.

## References

[1]: M. Soltani, C. A. Vargas-Garcia, and A. Singh, "Conditional Moment Closure Schemes for Studying Stochastic Dynamics of Genetic Circuits", IEEE Transactions on Biomedical Circuits and Systems 9, 518–526 (2015). [https://doi.org/10.1109/tbcas.2015.2453158](https://doi.org/10.1109/tbcas.2015.2453158)

[2]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Approximation and inference methods for stochastic biochemical kinetics - a tutorial review", Journal of Physics A: Mathematical and Theoretical 50, 093001 (2017). [https://doi.org/10.1088/1751-8121/aa54d9](https://doi.org/10.1088/1751-8121/aa54d9)