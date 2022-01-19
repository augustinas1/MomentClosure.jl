# [Moment Expansion (SDE)](@id moment_expansion_SDE)

Moment equations can be also written down for models described by $N$-dimensional Stochastic Differential Equations (SDEs) of the form 
```math
\begin{align*}
    d\mathbf{n} = \mathbf{f} \left(\mathbf{n}, t \right) dt + \mathbf{G}\left( \mathbf{n}, t \right) d\mathbf{w}_t,
\end{align*}
```
where $\mathbf{n}=(n_1, \dotsc, n_N)$ is a $N \times 1$ state vector, $\mathbf{f} \left(\mathbf{n}, t \right) = \left( f_1(\mathbf{n}, t), \dotsc, f_N(\mathbf{n}, t) \right)^\top$ is the deterministic (drift) term, $\mathbf{G}(\mathbf{n}, t)$ is a $N \times m$ (diffusion) matrix and $\mathbf{w}_t$ is a $m \times 1$ vector of [Wiener processes](https://en.wikipedia.org/wiki/Wiener_process) so that $\langle d\mathbf{w}_t \rangle = \mathbf{0}$ and $d\mathbf{w}_t d\mathbf{w}_t^\top = \mathbf{I}\,dt$, with $\mathbf{I}$ being an $m \times m$ identity matrix. Note that the notation used throughout is consistent with that of the [previous section](@ref moment_expansion_CME), and the state vector here is the same as for a chemical reaction network of $N$ species (the only difference being that its elements are now continuous variables). The derivations below are based on Ghusinga et al. (2017) [1] and Bover (1978) [2].

## [Raw Moment Equations](@id raw_moment_eqs_SDE)

To obtain the raw moment equations, we begin with [Itô's lemma](https://en.wikipedia.org/wiki/It%C3%B4%27s_lemma) stating that for any smooth scalar function $h(\mathbf{n})$ one has 
```math
\begin{align*}
    dh(\mathbf{n}) = \left\{ \sum_{i=1}^N \frac{\partial h(\mathbf{n}) }{\partial n_i} f_i + \frac{1}{2} \sum_{i=1}^N \sum_{j=1}^N \frac{\partial^2 h(\mathbf{n}) }{\partial n_i \partial n_j} (\mathbf{G}\mathbf{G}^\top)_{ij} \right\}dt + \sum_{i=1}^N \sum_{j=1}^N \frac{\partial h(\mathbf{n}) }{\partial n_i} \mathbf{G}_{ij} dw_j \,.
\end{align*}
```
Taking the expectation of both sides we obtain
```math
\begin{align*}
    \frac{d \langle h(\mathbf{n}) \rangle}{dt} = \sum_{i=1}^N \left\langle \frac{\partial h(\mathbf{n})}{\partial n_i} f_i \right\rangle + \frac{1}{2} \sum_{i=1}^N \sum_{j=1}^N \left\langle  \frac{\partial^2 h(\mathbf{n})}{\partial n_i \partial n_j} (\mathbf{G}\mathbf{G}^\top)_{ij}\right\rangle \,,
\end{align*}
```
where we have used the fact that $\langle d\mathbf{w} \rangle = \mathbf{0}$. Finally, the general form of raw moment equations up to $m^{\text{th}}$ order is obtained by setting $h(\mathbf{n}) = \mathbf{n}^{\mathbf{i}} = n_1^{i_1} \dotsm n_N^{i_N}$, so that $\langle h(\mathbf{n}) \rangle = \mu_{\mathbf{i}}$ and $|\mathbf{i}| = \sum_{j=1}^N i_j \leq m$, leading to 
```math
\begin{align*}
    \frac{d\mu_{\mathbf{i}}}{dt} = \sum_{j=1}^N \left\langle \frac{\partial \mathbf{n}^{\mathbf{i}} }{\partial n_j} f_j \right\rangle + \frac{1}{2} \sum_{j=1}^N \sum_{k=1}^N \left\langle \frac{\partial^2 \mathbf{n}^{\mathbf{i}} }{\partial n_j \partial n_k} (\mathbf{G}\mathbf{G}^\top)_{jk}\right\rangle \,.
\end{align*}
```
The equation for the *mean* vector $\mathbf{μ} = \langle \mathbf{n} \rangle = (\mu_1, \dotsc, \mu_N)$ take a particularly simple form (as the second derivatives are zero):
```math
\begin{align*}
    \frac{d\mathbf{μ}}{dt} = \langle \mathbf{f}(\mathbf{n}, t) \rangle.
\end{align*}
```

## [Central Moment Equations](@id central_moment_eqs_SDE)

Note that the raw moment expansion formulated above is valid only when the SDE terms $\mathbf{f} \left(\mathbf{n}, t \right)$ and $\mathbf{G}(\mathbf{n}, t)$ are *polynomial* functions. If these functions are *non-polynomial*, we have to consider central moment equations and Taylor expand both $\mathbf{f} \left(\mathbf{n}, t \right)$ and $\mathbf{G}(\mathbf{n}, t)$ around the mean $\mathbf{μ}$ up to the expansion order $q$ ([as done in case of the CME with non-polynomial propensities](@ref central_moment_eqs)).

The corresponding Taylor expansions can be written down as:
```math
\begin{align*}
    f_j(\mathbf{n}) &= \sum_{|\mathbf{k}|=0}^{q} \frac{1}{\mathbf{k}!}D^{\mathbf{k}} f_j(\mathbf{μ}) (\mathbf{n} - \mathbf{μ})^{\mathbf{k}} + \dotsb \,, \\
    \left( \mathbf{G}(\mathbf{n}) \mathbf{G}(\mathbf{n})^\top \right)_{jk} &= \sum_{|\mathbf{l}|=0}^{q} \frac{1}{\mathbf{l}!} D^{\mathbf{l}} \left( \mathbf{G}(\mathbf{μ})\mathbf{G}(\mathbf{μ})^\top \right)_{jk} (\mathbf{n} - \mathbf{μ})^{\mathbf{l}} + \dotsb \,,
\end{align*}
```
and we remind the reader that we have [previously](@ref central_moment_eqs) introduced the notation
```math
\begin{align*}
\mathbf{k}! &= k_1!\dotsm k_N!, \\
D^{\mathbf{k}} f &= \frac{\partial^{|\mathbf{k}|}}{\partial n_1^{k_1} \dotsm \partial n_N^{k_N}} f,\\
(\mathbf{n} - \mathbf{μ})^{\mathbf{k}} &= (n_1-\mu_1)^{k_1} \dotsm (n_N-\mu_N)^{k_N} \;.
\end{align*}
```
It follows that the equations for the means can now be written down as:
```math
\begin{align*}
    \frac{d\mu_i}{dt} = \langle f_i \rangle = \sum_{|\mathbf{j}|=0}^q \frac{1}{\mathbf{j}!} D^{\,\mathbf{j}} f_i(\mathbf{μ}) M_{\mathbf{j}} + \dotsb \,,
\end{align*}
```
where $M_{\mathbf{j}} = \langle (\mathbf{n}-\mathbf{μ})^{\mathbf{j}} \rangle$ are central moments.

To obtain the central moment equations, we first perform a change of variables $\mathbf{y} = \mathbf{n} - \langle \mathbf{n} \rangle$ and note again that $d\langle \mathbf{n} \rangle/dt = \langle \mathbf{f}(\mathbf{n}, t) \rangle$, so that the SDE becomes
```math
\begin{align*}
    d\mathbf{y} = \left[ \mathbf{f}\left(\mathbf{y} + \langle\mathbf{n}\rangle, t \right) - \left\langle \mathbf{f} \left(\mathbf{y} + \langle\mathbf{n}\rangle, t \right) \right\rangle \right] dt + \mathbf{G}\left(\mathbf{y} + \langle\mathbf{n}\rangle, t \right) d\mathbf{w}_t \,.
\end{align*}
```
Then, using Itô's lemma for a smooth function $h(\mathbf{y})$ and taking the expectation of both sides we find
```math
\begin{align*}
    \frac{d \langle h(\mathbf{y}) \rangle}{dt} = \sum_{i=1}^N \left[ \left\langle \frac{\partial h(\mathbf{y})}{\partial y_i} f_i \right\rangle - \left\langle \frac{\partial h(\mathbf{y})}{\partial y_i} \right\rangle \left\langle f_i \right\rangle \right] + \frac{1}{2} \sum_{i=1}^N \sum_{j=1}^N \left\langle  \frac{\partial^2 h(\mathbf{y})}{\partial y_i \partial y_j} (\mathbf{G}\mathbf{G}^\top)_{ij}\right\rangle \,.
\end{align*}
```
Letting $h(\mathbf{y}) = (\mathbf{n} - \mathbf{μ})^{\mathbf{i}}$ and substituting in the needed Taylor expansions we obtain the final form of central moment equations:
```math
\begin{align*}
    \frac{d M_{\mathbf{i}}}{dt} &= \sum_{j=1}^N \Bigg[ \left\langle \frac{\partial (\mathbf{n} - \mathbf{μ})^{\mathbf{i}} }{\partial (n_j - \mu_j)} \sum_{|\mathbf{k}|=0}^{q-|\mathbf{i}|+1} \frac{1}{\mathbf{k}!} D^{\mathbf{k}} f_j(\mathbf{μ}) (\mathbf{n}-\mathbf{μ})^\mathbf{k} + \dotsb \right\rangle \\
    &- \left\langle \frac{\partial (\mathbf{n}-\mathbf{μ})^{\mathbf{i}} }{\partial(n_j - μ_j)} \right\rangle \left\langle \sum_{|\mathbf{k}|=0}^{q} \frac{1}{\mathbf{k}!} D^{\mathbf{k}} f_j(\mathbf{μ}) (\mathbf{n}-\mathbf{μ})^\mathbf{k} + \dotsb \right\rangle \Bigg] \\
    &+ \frac{1}{2} \sum_{j=1}^N \sum_{k=1}^N \left\langle \frac{\partial^2 (\mathbf{n}-\mathbf{μ})^{\mathbf{i}}}{\partial(n_j-\mu_j)\partial(n_k-\mu_k)} \sum_{|\mathbf{l}|=0}^{q-|\mathbf{i}|+2} \frac{1}{\mathbf{l}!} D^{\,\mathbf{l}} \left( \mathbf{G}(\mathbf{μ}) \mathbf{G}(\mathbf{μ})^\top \right)_{jk} (\mathbf{n}-\mathbf{μ})^\mathbf{l} + \dotsc \right\rangle
\end{align*}
```

## References

[1]: K. R. Ghusinga, M. Soltani, A. Lamperski, S. V. Dhople, and A. Singh, "Approximate moment dynamics for polynomial and trigonometric stochastic systems", IEEE 56th Annual Conference on Decision and Control (2017). [https://doi.org/10.1109/CDC.2017.8263922](https://doi.org/10.1109/CDC.2017.8263922)

[2]: D. C. C. Bover, "Moment Equation Methods for Nonlinear Stochastic Systems", Journal of Mathematical Analysis and Applications 65, 306-320 (1978). [https://doi.org/10.1016/0022-247X(78)90182-8](https://doi.org/10.1016/0022-247X(78)90182-8)