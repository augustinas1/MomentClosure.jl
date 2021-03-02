# [Moment Expansion](@id moment_expansion)

## The Chemical Master Equation

Consider a chemical reaction network with $N$ different molecular species $X_i$ ($i=1, \dotsc, N$) and $R$ chemical reactions, so that the system can be described by [1]:

```math
\begin{align*}
    \sum_{i=1}^{N} s_{ij} X_{i} \xrightarrow{k_j} \sum_{i=1}^N r_{ij} X_i, \quad j=1,\dotsc,R,
\end{align*}
```
where $s_{ij}$ and $r_{ij}$ respectively denote the numbers of reactant and product molecules of species $i$ in the chemical reaction $j$. The stoichiometric matrix is defined as $S_{ij} = r_{ij} - s_{ij}$ (the net change in the number of molecules of species $X_{i}$ when the $r^{\text{th}}$ reaction occurs). The state of the system is determined by the state vector $\mathbf{n}=(n_1, \dotsc, n_N)$, where $n_i$ is the number of $X_i$ molecules. The time evolution of the probability distribution of $\mathbf{n}$ is described by the Chemical Master Equation (CME)
```math
\begin{equation}
    \frac{d P(\mathbf{n}, t)}{dt} = \sum^{R}_{r=1} \Big[ a_r(\mathbf{n}-S_r)P(\mathbf{n}-S_r, t) - a_r(\mathbf{n}) P(\mathbf{n}, t) \Big]\;, \tag{1}
\end{equation}
```
where $a_{r}(\mathbf{n})$ is the propensity function of the $r^{\text{th}}$ reaction and $S_r$ is the $r^{\text{th}}$ column of the stoichiometric matrix $S$.

Modelling using the CME framework is common in the study of biochemical and gene networks within cells. Although such master equations are rather simple in structure, for most systems there are no known analytical solutions and their stochastic simulations can be very computationally expensive. One possible approach to investigate the system at hand is to approximate the whole probability distribution solution of the CME in terms of its first few moments (e.g. mean and variance).

## [Raw Moment Equations](@id raw_moment_eqs)
From Eq. ([1](#mjx-eqn-1)) we can obtain a system of ordinary differential equations (ODEs) governing the time-evolution of the *raw* moments of the system up to specified order $m$, given by $\mu_{\mathbf{i}} = \mu_{i_1, \dotsc, i_N} = \langle n_1^{i_1} \dotsm n_N^{i_N} \rangle$, where the indices $\mathbf{i}=(i_1, \dotsc, i_N)$ are restricted such that $|\mathbf{i}| = \sum_{j=1}^N i_j \leq m$. Note that first order raw moments are simply the means. For example, when N = 3, the mean molecule number of second species in the system is given by $\mu_2 = \mu_{0, 1, 0} = \langle n_2 \rangle$ (when $N=3$). We chose to relax the notation throughout so that the means can be indicated by a single index of the molecular species (used in this section) or by the corresponding unit vector (used in the code).

To find the first moment equations, we multiply the CME by $n_i$ and sum over all possible states:
``` math
\begin{align*}
    \sum_{\mathbf{n}} n_i \frac{dP(\mathbf{n},t)}{dt} &= \sum_{n_1}^{\infty} \sum_{n_2}^{\infty} \dotsm \sum_{n_N}^{\infty} n_i \frac{dP(\mathbf{n},t)}{dt} \\
    &= \sum_{r} \sum_{\mathbf{n}} n_i a_r(\mathbf{n}-S_r)P(\mathbf{n}-S_r, t) - n_i a_r(\mathbf{n}) P(\mathbf{n}, t) \;.
\end{align*}
```
Applying a transformation $\mathbf{n}-S_r \rightarrow \mathbf{n}$ on the first term in the sum leads to
``` math
\begin{align*}
    \frac{d\mu_i}{dt} &= \sum_{r} \sum_{\mathbf{n}} (n_i+S_{ir}) a_r(\mathbf{n})P(\mathbf{n}, t) - n_i a_r(\mathbf{n}) P(\mathbf{n}, t) \\
    &= \sum_{r} \sum_{\mathbf{n}} S_{ir} a_r({\mathbf{n}}) P(\mathbf{n}, t) \\
    &= \sum_{r} S_{ir} \langle a_r(\mathbf{n}) \rangle\;. \tag{2}
\end{align*}
```
The derivation can be extended to the multivariate case of higher order moments $\mu_{\mathbf{i}}$ of arbitrary order $m$:
``` math
\begin{align*}
    \frac{\mu_{\mathbf{i}}}{dt} &= \sum_{\mathbf{n}} n_1^{i_1}\dotsm n_N^{i_N} \frac{dP(\mathbf{n}, t)}{dt} \\
    &= \sum_{r} \sum_{\mathbf{n}} \left[ (n_1+S_{1r})^{i_1}\dotsm (n_N +S_{2r})^{i_N} - n_1^{i_1}\dotsm n_N^{i_N} \right] a_r(\mathbf{n})P(\mathbf{n}, t) \\
    &= \sum_r \sum_{\mathbf{n}} \left[ \sum^{i_1}_{j_1=0} \binom{i_1}{j_1} n_{1}^{j_1}S_{1r}^{\,i_1-j_1} \dotsm \sum^{i_N}_{j_N=0} \binom{i_N}{j_N} n_{N}^{j_N}S_{Nr}^{\, i_N-j_N} - n_1^{i_1}\dotsm n_N^{i_N} \right] a_r(\mathbf{n})P(\mathbf{n}, t) \\
    &= \sum_r \sum_{\mathbf{n}} \sum^{|\mathbf{i}|-1}_{|\mathbf{j}|=0} \binom{i_1}{j_1} \dotsm \binom{i_N}{j_N} S_{1r}^{\,i_1-j_1} \dotsm S_{Nr}^{\, i_N-j_N} n_{1}^{j_1} \dotsm n_{N}^{j_N} a_r(\mathbf{n})P(\mathbf{n}, t) \\
    &= \sum_r \sum^{|\mathbf{i}|-1}_{|\mathbf{j}|=0} \binom{\, \mathbf{i} \,}{\, \mathbf{j} \,} S_r^{\,\mathbf{i}-\mathbf{j}} \langle \mathbf{n}^{\mathbf{j}} a_r(\mathbf{n}) \rangle \;,
\end{align*}
```
where we have introduced multi-index notation:
```math
\begin{align}
    |\mathbf{i}| &= \sum^{N}_{j=1} i_j , \tag{3} \\
     \sum^{|\mathbf{i}|-1}_{|\mathbf{j}|=0} &= \sum_{\substack{0 \leq j_1 \leq i_1, \, \dotsc, \, 0 \leq j_N \leq i_N\\ 0 \leq j_1+\dotsb+j_N \leq |\mathbf{i}|-1}} \tag{4} \\
    \binom{\, \mathbf{i} \,}{\, \mathbf{j} \,} &= \binom{i_1}{j_1} \dotsm \binom{i_N}{j_N} \tag{5}, \\
    \mathbf{n}^{\mathbf{j}} &= n_1^{j_1}\dotsm n_N^{j_N}, \tag{6} \\
    S_{r}^{\,\mathbf{i}-\mathbf{j}} &= S_{1r}^{\, i_1-j_1} \dotsm S_{Nr}^{\, i_N-j_N} \tag{7}\;.
\end{align}
```
**It is crucial to stress three key points:**
- Throughout all derivations presented in this section we assume that the components of the net stoichiometry matrix $S$ are constant values. However, in certain cases the reaction product may itself be a stochastic variable so that the corresponding expectation values, $\langle S_{r}^{\,\mathbf{i}-\mathbf{j}} \rangle$, must be taken into account. We have implemented a limited support for such systems where components of $S$ can be *independent geometrically distributed* variables, see the specific [gene network example](@ref gene_network_example).
- Generation of raw moment equations is only possible if the kinetics of the system at hand are governed by the law of mass action where *all* propensity functions are *polynomials* in $\mathbf{n}$, otherwise the expectation value terms $\langle a_r(\mathbf{n}) \rangle$ and $\langle n_{i_1}^{j_1} \dotsm n_{i_N}^{j_N} a_r(\mathbf{n}) \rangle$ are ill-defined. This key issue can be overcome by the general central moment expansion method presented below.
- The order of the polynomials $a_r(\mathbf{n})$ determines the order of moments encountered in the generated system of ODEs. If the system is linear (contains only zeroth and first order reactions), the $m^{\text{th}}$ order moments will depend only of moments of order $m$ or lower, hence constituting a *finite* hierarchy of moment ODEs that can be readily solved numerically or otherwise without approximations. However, if the system is non-linear (involves second or higher order reactions), moment equations will depend on higher order moments. For example, if the reaction network contains bimolecular reactions, the corresponding propensity functions will be second order polynomials and hence $m^{\text{th}}$ order moment equations will now depend on $(m+1)^{th}$ order moments. This leads to an infinite hierarchy of coupled moment equations where each moment will depend on higher order moments—it cannot be solved directly and needs to be truncated. This can be achieved using one of many *moment closure approximations* that express $(m+1)^{\text{th}}$ order moments in terms of $m^{\text{th}}$ and lower order moments using different distributional assumptions, effectively closing the hierarchy and enabling one to solve the moment equations up to $m^{\text{th}}$ order [1]. Details of all closure methods currently implemented within the package can be found in the next **Theory** [section](@ref moment_closure_approximations).

## [Central Moment Equations](@id central_moment_eqs)

As we have seen, the moment equations of $P(\mathbf{n}, t)$ can be obtained in a straightforward manner if the kinetics of the system are governed by the law of mass action where all propensity functions are polynomials in $\mathbf{n}$ [1]. Similarly, given all propensities are rational functions, a polynomial form can also be recovered [2]. However, problems arise when the propensities take more complicated non-polynomial functions. This can nevertheless be overcome by considering a more general method of moment expansion that enables us to obtain mean and central moment equations up to arbitrary order for virtually any chemical reaction network with any type of *smooth* (infinitely differentiable) propensity functions. Such framework was first independently formulated by Lee [3] and Ale et al. [4]—our derivation below closely follows these works.

We start by Taylor-expanding the propensity functions around the mean $\mathbf{μ} = \langle \mathbf{n} \rangle = (\mu_1, \dotsc, \mu_N)$. This allows us to consider any general propensity function under the assumption that it is infinitely differentiable (smooth). The expansion leads to
```math
\begin{align*}
    a_r({\mathbf{n}}) &= a_{r}(\mathbf{μ}) + \sum_{j_1}^N \frac{\partial a_r(\mathbf{μ})}{\partial n_{j_1}}(n_{j_1}-\mu_{j_1}) \\
    &+ \frac{1}{2!}\sum_{j_1, j_2} \frac{\partial^2 a_r({\mathbf{μ}})}{\partial n_{j_1} \partial n_{j_2}} (n_{j_1}-\mu_{j_1}) (n_{j_2}-\mu_{j_2}) + \dotsb \\
    & + \frac{1}{q!} \sum_{j_1, j_2, \dotsc, j_q} \frac{\partial^q a_r(\mathbf{μ})}{\partial n_{j_1} \partial n_{j_2}\dotsm\partial n_{j_q}} (n_{j_1}-\mu_{j_1}) (n_{j_2}-\mu_{j_2}) \dotsm(n_{j_q} - \mu_{j_q}) + \dotsb \\
    &= \sum_{|\mathbf{j}|=0}^{q} \frac{1}{\mathbf{j}!} D^{\mathbf{j}} a_r(\mathbf{μ}) (\mathbf{n} - \mathbf{μ})^{\mathbf{j}} + \dotsb \;,
\end{align*}
```
where $q$ controls the expansion order and we have simplified the expression by using Eqs. ([3-7](#mjx-eqn-3)) in addition to:
```math
\begin{align*}
\mathbf{j}! &= j_1!\dotsm j_N!, \\
D^{\,\mathbf{j}} f &= \frac{\partial^{|\mathbf{j}|}}{\partial n_1^{j_1} \dotsm \partial n_N^{j_N}} f \;.
\end{align*}
```
Now we can obtain equations governing the time-evolution of the means $μ_i$ and the central moments of arbitrary order $m$, $M_{\mathbf{i}} = M_{i_1, \dotsc, i_N} = \langle (n_1 - \mu_1)^{i_1}\dotsm(n_N - \mu_N)^{i_N} \rangle$, where again $|\mathbf{i}| = m$. We first consider the equations for the means immediately starting from Eq. ([2](#mjx-eqn-2)):
``` math
\begin{align*}
    \frac{d\mu_i}{dt} &= \sum_{r} S_{ir} \langle a_r(\mathbf{n}) \rangle \\
    &= \sum_r S_{ir} \Big( \sum_{|\mathbf{j}|=0}^{q} \frac{1}{\mathbf{j}!} D^{\, \mathbf{j}} a_r(\mathbf{μ}) \langle (\mathbf{n} - \mathbf{μ})^{\mathbf{j}} \rangle + \dotsb \Big) \\
    &= \sum_r S_{ir} \Big( \sum_{|\mathbf{j}|=0}^{q} \frac{1}{\mathbf{j}!} D^{\, \mathbf{j}} a_r(\mathbf{μ}) M_{\mathbf{j}} + \dotsb \Big) \;.
\end{align*}
```
We can now derive the equations for central moments by taking the time derivative of $M_{\mathbf{i}}$ and using Eq. ([1](#mjx-eqn-1)):
```math
\begin{align*}
    \frac{dM_{\mathbf{i}}}{dt} &= \sum_{\mathbf{n}} (n_1 - \mu_1)^{i_1}\dotsm(n_N - \mu_N)^{i_N} \frac{dP(\mathbf{n}, t)}{dt} \\
    &- \sum_{\mathbf{n}} \sum_{j=1}^{N} i_j \frac{d\mu_j}{dt} (n_1-\mu_1)^{i_1} \dotsm(n_j-\mu_j)^{i_j-1}\dotsm(n_N-\mu_N)^{i_N}P(\mathbf{n}, t) \\
    &= \sum_{r} \sum_{\mathbf{n}} \Big[ (n_1 - \mu_1 + S_{1r})^{i_1} \dotsm (n_N-\mu_N+S_{Nr})^{i_N} \\
    &- (n_1-\mu_1)^{i_1}\dotsm(n_N-\mu_N)^{i_N} \Big]a_r(\mathbf{n})P(\mathbf{n}, t) - \sum_j i_j \frac{d\mu_{j}}{dt}M_{i_1,\dotsc,i_{j}-1,\dotsc, i_N} \\
    &= \sum_r \sum_{\mathbf{n}} \Bigg[ \sum_{|\mathbf{j}|=0}^{|\mathbf{i}|-1} \binom{i_1}{j_1} \dotsm \binom{i_N}{j_N} S_{1r}^{\, i_1-j_1} \dotsm S_{Nr}^{\, i_N-j_N}(n_1-\mu_1)^{j_1}\dotsm(n_N-\mu_N)^{j_N}\Bigg] \\
    &\times \Bigg[ \sum_{|\mathbf{j}|=0}^{q} \frac{1}{\mathbf{j}!} D^{\mathbf{j}} a_r(\mathbf{μ}) (\mathbf{n} - \mathbf{μ})^{\mathbf{j}} + \dotsb \Bigg]-\sum_j i_j \frac{d\mu_{j}}{dt} M_{\mathbf{i}_{j-}} \\
    &= \sum_r \sum_{|\mathbf{j}|=0}^{|\mathbf{i}|-1} \Bigg[ \binom{i_1}{j_1} \dotsm \binom{i_N}{j_N} S_{1r}^{\, i_1-j_1} \dotsm S_{Nr}^{\, i_N-j_N} \sum_{|\mathbf{k}|=0}^{q-|\mathbf{j}|} \frac{1}{\mathbf{k}!}D^{\mathbf{k}} a_r(\mathbf{μ}) M_{\mathbf{j}+\mathbf{k}} + \dotsb \Bigg] \\  &- \sum_j i_j \frac{d\mu_{j}}{dt} M_{\mathbf{i}_{j-}} \\
    &= \sum_r \sum_{|\mathbf{j}|=0}^{|\mathbf{i}|-1} \Bigg[ \binom{\mathbf{i}}{\mathbf{j}} S_r^{\mathbf{i}-\mathbf{j}} \sum_{|\mathbf{k}|=0}^{q-|\mathbf{j}|} \frac{1}{\mathbf{k}!}D^{\mathbf{k}} a_r(\mathbf{μ}) M_{\mathbf{j}+\mathbf{k}} + \dotsb \Bigg] - \sum_j i_j \frac{d\mu_{j}}{dt} M_{\mathbf{i}_{j-}} \;,
\end{align*}
```
where we have also defined
```math
\begin{align*}
    M_{\mathbf{i}_{j-}}&=M_{i_1,\dotsc,i_{j}-1,\dotsc, i_N}, \\
    M_{\mathbf{j}+\mathbf{k}} &= M_{j_1+k_1,\dotsc,j_N+k_N} \;.
\end{align*}
```
Although raw moment equations for non-linear mass-action systems already require approximate treatment using moment closure, here we have an additional complication: if a system contains *non-polynomial* propensity functions, the equations for both means and central moments will in principle depend on an *infinite* number of higher order central moments. Hence the *expansion order* $q$ is of utmost importance as it controls the degree of approximation: the Taylor expansion of propensity functions is performed up to $q^{\text{th}}$ order so that $m^{\text{th}}$ central moment equations will depend on central moments of order $q>m$ and lower. Finally, moment closure approximations can be applied similarly as in the case of raw moment equations.


## References

[1]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Approximation and inference methods for stochastic biochemical kinetics - a tutorial review", Journal of Physics A: Mathematical and Theoretical 50, 093001 (2017). [https://doi.org/10.1088/1751-8121/aa54d9](https://doi.org/10.1088/1751-8121/aa54d9)

[2]: P. Milner, C. S. Gillespie, and D. J. Wilkinson, "Moment closure approximations for stochastic kinetic models with rational rate laws", Mathematical Biosciences 231, 99-104 (2011). [https://doi.org/10.1016/j.mbs.2011.02.006](https://doi.org/10.1016/j.mbs.2011.02.006)

[3]: A. Ale, P. Kirk, and M. P. H. Stumpf, "A general moment expansion method for stochastic kinetic models", The Journal of Chemical Physics 138, 174101 (2013). [https://doi.org/10.1063/1.4802475](https://doi.org/10.1063/1.4802475)

[4]: C. H. Lee, "A Moment Closure Method for Stochastic Chemical Reaction Networks with General Kinetics", MATCH Communications in Mathematical and in Computer Chemistry 70, 785-800 (2013). [https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_785-800.pdf](https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_785-800.pdf)
