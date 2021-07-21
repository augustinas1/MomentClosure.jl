# [Linear Mapping Approximation](@id linear_mapping_approximation)

The Linear Mapping Approximation (LMA) provides a novel way of approximating the solution of the [CME](@ref chemical_master_equation) and has been shown to be accurate for a variety of models of gene regulatory networks (GRNs) [1]. It is based on mapping a *nonlinear* GRN onto an equivalent *linear* GRN so that the exact solution of the linear system gives an approximate solution of the nonlinear system. The LMA is restricted in its applicability to chemical reaction networks where one of the substrates in each nonlinear reaction is a molecular species which copy number can be either zero or one (a Bernoulli/binary random variable). Note that a network can contain an arbitrary number of such species but more than one of them cannot be involved in any nonlinear reaction. Below we provide a short overview of the LMA and urge the reader to see the original [paper](https://doi.org/10.1038/s41467-018-05822-0) for a more comprehensive description [1].

To explain the LMA procedure, we start by considering a simple example of a two-state GRN as done in a [previous section on conditional closures](@ref conditional_closures), denoting the binary gene state by $g$ and the protein number by $p$. We assume that there is a single nonlinear reversible reaction in the network given by
```math
\begin{align*}
G+P &\underset{σ_u}{\stackrel{σ_b}{\rightleftharpoons}} G^*,
\end{align*}
```
where $P$ denotes the protein and the gene can be in either state $G$ ($g=1$) or $G^*$ ($g=0$). Our aim is to find an approximate time-dependent probability distribution of protein numbers $p$ at time $t$.

The steps of performing the LMA can then be described as follows:

**1.** Find the linear network by replacing any reversible nonlinear reaction (it must involve one binary species) in the nonlinear network by a reversible pseudo first-order reaction between the binary species' states.
- In our example, we replace the reaction above with $G \underset{σ_u}{\stackrel{\bar{σ}_b}{\rightleftharpoons}} G^*$, noting that the rate parameter is changed from $σ_b$ to $\bar{σ}_b$. Moreover, this approach is just as applicable in case of cooperativity, e.g., $G+nP \rightleftharpoons G^*$ (where $n$ is an integer indicating the cooperative order) would be similarly replaced with $G \rightleftharpoons G^*$.

**2.** Approximate the changed rate parameters of the linearised reactions by their expectation values.
- As noted in [1], the first-order reaction $G \stackrel{\bar{σ}_b}\rightarrow G^*$ maps onto the second-order reaction $G+P \stackrel{σ_b}\rightarrow  G^*$ if we choose $\bar{σ}_b = σ_b \left(p \,|\, g=1 \right)$, where $p \,|\, g=1$ indicates the instantaneous protein number given the gene is in the state $G$. In LMA, we use the mean-field approximation taking the expectation value of the rate so that
  ```math
  \begin{align*}
  \bar{σ}_b = σ_b \left \langle p \,|\, g=1 \right \rangle = σ_b \frac{\left\langle pg \right\rangle}{\left\langle g \right\rangle}.
  \end{align*}
  ```
  The same procedure can be extended to the general nonlinear reaction where $n$ proteins bind cooperatively. The effective parameter is then given by
  ```math
  \bar{σ}_b = σ_b \frac{\left\langle \prod_{i=0}^{n-1} \left( p-i\right)g \right\rangle}{\left\langle g \right\rangle}.
  ```

**3.** Write down the moment equations for the linear network using the approximated stochastic rates.
- Note that the moment equations must be generated up the order given by the highest order nonlinear reaction in the network. If the only nonlinear reaction is the second-order reaction $G+P \stackrel{σ_b}\rightarrow G^*$, we need to consider only moments up to the second order (as hinted by the functional form of $\bar{σ}_b$ above). However, the moment hierarchy is otherwise closed, no additional moment closure approximations need to be performed, and therefore we can solve the moment equations in a straightforward manner.

**4.** Solve the moment equations numerically up to time $t$ and plug the resulting moment values into the equations for the effective parameters. Proceed to calculate the time-average of these parameters over the time-interval $[0, t]$.
- In our example, plugging the solved-for moment values into the equation for $\bar{σ}_b$ allows us to interpret the effective parameter as a time-dependent function $\bar{σ}_b(t)$. However, as the time-dependent probability distribution solution of the CME for the nonlinear network with a general time-dependent $\bar{σ}_b$ is most likely intractable, we transform $\bar{σ}_b(t)$ into a time-independent constant by taking its time-average $\bar{σ}_b^* = \int_0^t \bar{σ}_b(t') dt' / t$. This approach is justified in [1] by considering the Magnus expansion.

**5.** Obtain the time-dependent probability distribution solution of the CME of the linear network assuming that the rate parameters of the linearised reactions are time-independent constants.
- Note that this step is the major limitation of the LMA as closed-form solutions are available only for a handful of systems (consult [1] for more details).

**6.** Finally, construct the approximate probability distribution of the nonlinear network at time $t$ by replacing the respective rate parameters with their time-averaged equivalents obtained in the previous step.

MomentClosure.jl provides automated generation of the closed moment equations using LMA given a nonlinear chemical reaction network and its linear equivalent. This encapsulates the first three steps of the LMA procedure outlined above which are general and can be seen as an original moment closure approximation. We apply the LMA on simple models of nonlinear GRNs and also discuss how the subsequent LMA steps can be performed in Julia on a case-by-case basis in [this tutorial example](@ref linear_mapping_approximation_example).

## References

[1]: Z. Cao and R. Grima, "Linear mapping approximation of gene regulatory networks with stochastic dynamics", Nature Communications 9, 3305 (2018). [https://doi.org/10.1038/s41467-018-05822-0](https://doi.org/10.1038/s41467-018-05822-0)
