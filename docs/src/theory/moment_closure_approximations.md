# Moment Closure Approximations

In the [previous section](moment_expansion.md), we have shown that for a non-linear system an infinite hierarchy of coupled moment equations is obtained that cannot be solved directly and therefore needs to be truncated in an approximate way. This can be achieved using moment closure approximations (MAs) using which all moments above a certain order $m$ are expressed in terms of $m^{\text{th}}$ and lower order moments using various (usually distributional) assumptions [1]. Doing so enables us to effectively *close* the moment hierarchy, leading to a finite set of ODEs which can then be solved numerically. In this section, we present some the commonly used MA methods that are implemented in MomentClosure. Please see the **Tutorials** section for examples showing different MAs applied to a variety of systems.

## Zero closure

The simplest MA is the "central-moment-neglect" MA (CMN-MA) [2], also referred to as "zero-closure" [3] or "low dispersion moment closure" [4], where CMN-MA at $m^{\text{th}}$ order means that the moment equations are truncated by setting all *central moments* above order $m$ to zero. For example, in the simple case of $2$nd order truncation, the moment equations for the means $\mu_i$ and covariances $C_{ij}$ become:
```math
\begin{align*}
    \frac{d \mu_i}{dt} &= \sum_{r} S_{ir} \Big( a_r(\mathbf{μ}) + \frac{1}{2}\sum_{i_1, i_2} \frac{\partial^2  a_r({\mathbf{μ}})}{\partial n_{i_1} \partial n_{i_2}} M_{i_1, i_2} \Big), \\
    \frac{dC_{ij}}{dt} &= \frac{dM_{\mathbf{0}_{i+, j+}}}{dt} \\
    &= \frac{d\langle (n_{i}-\mu_i) (n_{j}-\mu_{j}) \rangle}{dt} = \\
    &= \sum_r \Big[ S_{ir} \sum_k \frac{\partial a_r(\mathbf{μ})}{\partial n_k} C_{jk} + S_{jr} \sum_k \frac{\partial a_r(\mathbf{μ})}{\partial n_k} C_{ik} \\
    &+ S_{ir}S_{jr} \Big( a_r(\mathbf{μ}) + \frac{1}{2} \sum_{k, l} \frac{\partial^2 a_r(\mathbf{μ})}{\partial n_{k} \partial n_{l}} C_{kl} \Big) \Big] \;.
\end{align*}
```

## Normal closure

Another popular MA is the "normal moment-closure", pioneered by Goodman [5] and Whittle [6], where all *cumulants* $\kappa_{\mathbf{i}}$ above order $m$ are set to zero, approximating the probability distribution of the system with the normal distribution [2]:
```math
\begin{align*}
    \kappa_{\mathbf{i}} = 0, \quad \text{for} \; |\mathbf{i}| > m.
\end{align*}
```
In order to truncate the higher order central or raw moments $M_{\mathbf{i}}$ using normal closure we express them in terms of cumulants $\kappa_{\mathbf{i}}$ using the multivariate moment and cumulant relationships formalised by Balakrishnan et al. [7].

Note that a different implementation of normal closure can be found in literature [3], where the higher order central moments are expressed in terms of a sum of product of covariances using [Isserlis' theorem](https://en.wikipedia.org/wiki/Isserlis%27_theorem). However, one could argue that such formulation is not advisable as it assumes stronger "Gaussianity" of the underlying distribution than setting the higher order cumulants to zero which is less of an approximation on the form of the distribution and hence is preferable in the development of MAs. For example, in case we are truncating the moment equations at $4$th order, the truncation-order central moments would be expressed only in terms of covariances whereas our formulation using cumulants would explicitly include information about the computed values of third central moments, which is expected to improve numerical stability and lead to more accurate moment estimates.

## Poisson closure

Although the Poisson distribution lacks a general formulation for multiple variables [3], "Poisson MA" has been formulated [2, 8] assuming that the joint multivariate distribution is a product of univariate Poisson distributions, i.e., $n_i \sim \text{Poisson}(\mu_i)$. The cumulants of a univariate Poisson distribution are equal to the mean, hence in Poisson closure we set all higher order diagonal cumulants to the corresponding mean values and mixed cumulants to zero [2], which in our notation can be expressed as:
```math
\begin{align*}
    \kappa_{\mathbf{i}} &= \mu_k, \quad \text{if} \; |\mathbf{i}| > m \; \text{and} \; i_1,\dotsc,i_N = k \; \text{or} \; 0 \quad \text{for some} \; k \in \{1,\dotsc,N\}, \\
    \kappa_{\mathbf{i}} &= 0, \quad \; \; \text{if} \; |\mathbf{i}| > m \; \text{and} \; i_k \neq i_l \quad \text{for some} \; k,l \in \{1,\dotsc,N\}.
\end{align*}
```
Similarly to normal closure, the higher order central/raw moments can be expressed in terms of cumulants as described in [7].

## Log-normal closure

"Log-normal" MA, first applied by Keeling [9], allows one to truncate the moment equations under the assumption that the distribution of the underlying stochastic process is log-normal.
A positive multi-dimensional random variable $\mathbf{n}$ follows a log-normal distribution if its logarithm is normally distributed, so that $\mathbf{y} = \ln \mathbf{n}$, and $\mathbf{y} \sim \mathcal{N}(\mathbf{\nu}, \Sigma)$, where $\mathbf{\nu}$ and $\Sigma$ denote the vector of means and the covariance matrix respectively. By considering the moment generating function of the normal distribution, $\mathcal{N}(\mathbf{\nu}, \Sigma)$, one can show that the raw moments are given by [3, 10]:
```math
\begin{align*}
    \mu_{\mathbf{i}} = \exp\left(\mathbf{i}^{\top}\mathbf{\nu} + \frac{1}{2}\mathbf{i}^{\top}\Sigma\mathbf{i}\right).
\end{align*}
```
It follows that
```math
\begin{align*}
    \nu_i &= \ln \mu_i - \frac{1}{2} \Sigma_{ii}, \\
    \Sigma_{ij} &= \ln \left( 1 + \frac{\langle (n_i - \mu_i)(n_j - \mu_j) \rangle}{\exp \left( \nu_i + \nu_j + \frac{1}{2} \left( \Sigma_{ii} + \Sigma_{jj} \right) \right)} \right)  \\
    &= \ln \left( 1 + \frac{C_{ij}}{\exp \left( \nu_i + \nu_j + \frac{1}{2} \left( \Sigma_{ii} + \Sigma_{jj} \right) \right)} \right), \\
    \Sigma_{ii} &= \ln \left( 1 + \frac{C_{ii}}{\mu_i^2}  \right) .
\end{align*}
```
Note that central moments can be obtained from raw moments by utilising their general multivariate relationship [11].

## Gamma closure

The method of "gamma closure" was originally implemented by Lakatos et al. [3], where the authors acknowledged the ambiguity arising in defining multivariate gamma distributions, and, building upon previous definitions in the literature (e.g. [12] and [13]), proposed a new formulation of a multivariate gamma distribution. Here we reproduce the definition by closely following the description in [3] and elucidating some of the derivation steps.

We denote a random variable drawn from gamma distribution with shape $\alpha$ and scale $\beta$ as $n \sim \text{Gamma}(\alpha, \beta)$. The probability density function of the univariate gamma distribution with the corresponding shape-scale parameterisation is
```math
\begin{align*}
    f(n; \alpha, \beta) = \frac{1}{\Gamma(\alpha)\beta^{\alpha}} n^{\alpha-1}e^{-\frac{n}{\beta}},
\end{align*}
```
where $\Gamma$ is the gamma function. The $i^{\text{th}}$ raw moment of $n$ is given by
```math
\begin{align*}
    μ_i &= \frac{\Gamma(\alpha+i)\beta^i}{\Gamma(\alpha)} = (\alpha)_i \beta^i, \tag{a}
\end{align*}
```
where $(\alpha)_i = \alpha (\alpha+1) \dotsm (\alpha+i-1)$. Note that the moment generating function of X is
```math
\begin{align*}
    G_{n}(k) = \langle e^{kn} \rangle = \left( 1 - \beta k \right)^{-\alpha}.
\end{align*}
```
In order to construct a multivariate gamma distribution, we start by considering independent gamma variables $Y_{kl}$, $k, l = 1, \dotsc, N$, with shape and scale parameters $\alpha_{kl}$ and $\beta_{kl}$ respectively. Here we define $Y_{kl}$ to be symmetric in indices, i.e., $Y_{kl} = Y_{lk}$. Now consider an $N$-dimensional random vector $\mathbf{n} = \left( n_1, n_2, \dotsc, n_N \right)$, where $n_i$ is a linear combination of independent gamma variables:
```math
\begin{align*}
    n_i = \sum_{j=1}^N \frac{\beta_{ii}}{\beta_{ij}} Y_{ij}.
\end{align*}
```
The $i^{\text{th}}$ marginal moment generating function of the joint distribution of $\mathbf{n}$ is given by:
```math
\begin{align*}
    G_{n_i}(k_i) &= \left\langle e^{k_i \sum_{j=1}^N \frac{\beta_{ii}}{\beta_{ij}}} \right\rangle  \\
                 &= G_{Y_{i1}}\left(k_i \frac{\beta_{ii}}{\beta_{i1}} \right) G_{Y_{i2}}\left(k_i \frac{\beta_{ii}}{\beta_{i2}} \right) \dotsm G_{Y_{iN}}\left(k_i \frac{\beta_{ii}}{\beta_{iN}} \right)  \\
                 &= \left(1-\beta_{ii} k_i \right)^{\sum_{j=1}^N \alpha_{ij}},
\end{align*}
```
so that $n_i \sim \text{Gamma}(\alpha_i, \beta_i)$, where $\alpha_i = \sum_{j=1}^N \alpha_{ij}$ and $\beta_i = \beta_{ii}$. Therefore, we have obtained an $N$-variate gamma distribution, which can be denoted as $\mathbf{n} \sim MG(\mathbf{\alpha}, \mathbf{\beta})$, where the vectors of shape and scale parameters are given by $\mathbf{\alpha} = \left( \alpha_1, \dotsc, \alpha_N \right)$ and $\mathbf{\beta} = \left( \beta_1, \beta_2, \dotsc, \beta_N \right)$ respectively.

We can now readily obtain the raw $m^{\text{th}}$ order moment of $n_i$:
```math
\begin{align*}
    \langle n_i^{m} \rangle &= \left\langle \left( \sum_{j=1}^N \frac{\beta_{ii}}{\beta_{ij}} Y_{ij} \right)^m \right\rangle  \\
                            &= \left\langle \sum_{k_1+k_2+\dotsb+k_N = m} \frac{m!}{k_1!k_2! \dotsm k_N!} \prod_{j=1}^N \left( \frac{\beta_{ii}}{\beta_{ij}} Y_{ij} \right)^{k_j} \right\rangle  \\
                            &= \sum_{|\mathbf{k}|=m} \frac{m!}{\mathbf{k!}} \prod_{j=1}^N \left\langle \left( \frac{\beta_{ii}}{\beta_{ij}} Y_{ij} \right)^{k_j} \right\rangle  \\
                            &= \beta_i^m \sum_{|\mathbf{k}|=m} \frac{m!}{\mathbf{k!}} \left(
                            \prod_{j=1}^N \left( \alpha_{ij} \right)_{k_j} \right).
\end{align*}
```
Note that we have used Eq. ([a](#mjx-eqn-a)) to get to the last line. The mixed raw moments are computed in a similar fashion:
```math
\begin{align*}
    \mu_{\mathbf{i}} &= \left\langle \left( \sum_{j=1}^N \frac{\beta_{11}}{\beta_{1j}}Y_{1j} \right)^{i_1} \dotsm \left( \sum_{j=1}^N \frac{\beta_{NN}}{\beta_{Nj}}Y_{Nj} \right)^{i_N} \right\rangle  \\
                &= \mathbf{\beta}^{\mathbf{i}} \left\langle \sum_{|\mathbf{k_1}|=i_1} \dotsm \sum_{|\mathbf{k_N}|=i_N} \frac{\mathbf{i!}}{\mathbf{k_1!}\dotsm\mathbf{k_N!}} \prod_{j=1}^N \left(  \frac{Y_{1j}}{\beta_{1j}} \right)^{{k_1}_{j}} \dotsm \prod_{j=1}^N \left(  \frac{Y_{Nj}}{\beta_{Nj}} \right)^{{k_N}_{j}} \right\rangle  \\
                &= \mathbf{\beta}^{\mathbf{i}} \left\langle \sum_{|\mathbf{k_1}|=i_1} \dotsm \sum_{|\mathbf{k_N}|=i_N} \frac{\mathbf{i!}}{\mathbf{k_1!}\dotsm\mathbf{k_N!}} \prod_{q=1}^{N} \left( \frac{Y_{qq}}{\beta_{qq}} \right)^{{k_q}_{q}} \prod_{r=q+1}^N \left( \frac{Y_{qr}}{\beta_{qr}} \right)^{{k_{q}}_{r} + {k_r}_q} \right\rangle  \\
                 &= \mathbf{\beta}^{\mathbf{i}} \sum_{|\mathbf{k_1}|=i_1} \dotsm \sum_{|\mathbf{k_N}|=i_N} \frac{\mathbf{i!}}{\mathbf{k_1!}\dotsm\mathbf{k_N!}} \prod_{q=1}^{N} \left( \alpha_{qq} \right)_{{k_q}_{q}} \prod_{r=q+1}^N \left( \alpha_{qr} \right)_{{k_{q}}_{r} + {k_{r}}_{q}},
\end{align*}
```
where we have taken into account the symmetry in indices and defined each $\mathbf{k_i}$ as an $N$-dimensional vector $\mathbf{k_i} = \left( k_{i_{1}}, \dotsc, k_{i_{N}} \right)$ and $\mathbf{\beta}^{\mathbf{i}} = \beta_1^{i_1} \dotsm \beta_N^{i_N}$. Note that the mean and variance of $n_i$ can be obtained from:
```math
\begin{align*}
    \mu_i  &= \sum_{j=1}^N \alpha_{ij} \beta_i, \tag{b} \\
    C_{ii} &= (\alpha_i)_2 \beta_i^2 - \alpha_i^2\beta_i^2  \\
           &= \sum_{j=1}^N \alpha_{ij} \beta_i^2.
\end{align*}
```
Similarly, from Eq. ([a](#mjx-eqn-a)) it follows that
```math
\begin{align*}
    \langle n_i n_j \rangle = \beta_i \beta_j \left( \sum_{\substack{k, l \\ (k, l) \neq (j, i)}}^N  \alpha_{ik} \alpha_{jl} + (\alpha_{ij})_2 \right),
\end{align*}
```
which together with Eq. ([b](#mjx-eqn-b)) allows us to express the covariance as:
```math
\begin{align*}
    C_{ij} &= \langle n_i n_j \rangle - \beta_i \beta_j \left( \sum_{k,l} \alpha_{ik} \alpha_{jl} \right)  \\
           &= \alpha_{ij} \beta_i \beta_j.
\end{align*}
```
Finally, from the equations above we can obtain all shape and scale parameters:
```math
\begin{align*}
    \beta_i &= \frac{C_{ii}}{\mu_i}, \\
    \alpha_{ij} &= \frac{C_{ij}}{\beta_i \beta_j}, \\
    \alpha_{ii} &= \frac{\mu_i}{\beta_i} - \sum_{\substack{k \\ k \neq i}} \alpha_{ik}.
\end{align*}
```

## Derivative matching

The derivative matching MA [14, 15] is based on expressing moments above order $m$ in terms of lower order moments in such a way that their time derivatives match those of the *exact* moments at some initial time and initial conditions. We outline the approach below, closely following the complete exposition found in the original papers of Singh and Hespanha [14, 15].

The [raw moment equations](@ref raw_moment_eqs) up to order $m$ for any mass-action reaction network containing *at most* bimolecular (second order) reactions can be written down concisely in the matrix form:
```math
\begin{align*}
    \frac{d\mathbf{μ}}{dt} = \hat{\mathbf{a}} + A\mathbf{μ} + B \bar{\mathbf{μ}} \;,
\end{align*}
```
where $\mathbf{μ}$ is a vector containing all raw moments of the system up to order $m$ and $\bar{\mathbf{μ}}$ consists of all $(m+1)^{\text{th}}$ order raw moments which the equations depend on. The constant vector $\hat{\mathbf{a}}$ and constant matrices $A$ and $B$ are chosen appropriately for the system at hand. In this case, an MA can be defined as a procedure where each moment in $\bar{\mathbf{μ}}$, $\bar{μ}_{\mathbf{i}}$, is approximated by a *moment closure function* $\varphi_{\mathbf{i}} (\mathbf{μ})$ of moments up to order $m$. Then the moment equations can be rewritten as
```math
    \frac{d\mathbf{ν}}{dt} = \hat{\mathbf{a}} + A\mathbf{ν} + B\bar{\mathbf{φ}}(\mathbf{ν}) \;,
```
where the state of the system is now denoted by $\mathbf{ν}$ instead of $\mathbf{μ}$, stressing the fact that we are considering the approximation of the true moment dynamics, and $\bar{\mathbf{φ}}(\mathbf{ν})$ is the corresponding vector of moment closure functions.

The idea behind derivative matching is to determine a map $\bar{\mathbf{φ}}$ so that the time derivatives between the exact moments, $\mathbf{μ}(t)$, and the approximate moments, $\mathbf{ν}(t)$, would match at some initial time $t_0$ under the initial condition $\mathbf{μ}(t_0) = \mathbf{ν}(t_0)$:
```math
\begin{align*}
    \left. \frac{d^i \mathbf{μ}}{dt} \right\rvert_{t=t_0} = \left. \frac{d^i \mathbf{ν}}{dt} \right\rvert_{t=t_0}
\end{align*}
```
If these conditions hold, one can expect from a Taylor series approximation argument that $\mathbf{μ}(t)$ and $\mathbf{ν}(t)$ will stay close at least locally in time and hence the MA will be sufficiently accurate.

In order to move forward, Singh and Hespanha present what can be understood as essentially an ansatz. Firstly, moment closure functions for each $\mathbf{i}$, where $|\mathbf{i}| > m$, are chosen to have a separable form given by
```math
\begin{align*}
    φ_{\mathbf{i}}(\mathbf{μ}) = \prod_{\substack{1 \leq j_1+\dotsb+j_N \leq m}} \left( μ_{\mathbf{j}} \right)^{γ_{\mathbf{j}}} = \prod_{|\mathbf{j}|=1}^{m} \left( μ_{\mathbf{j}} \right)^{γ_{\mathbf{j}}} \;,
\end{align*}
```
where $γ_{\mathbf{j}}$ are constants (unique for each vector $\mathbf{i}$) that can be determined by solving the following linear equation system:
```math
\begin{align*}
    C^{\mathbf{i}}_{\mathbf{j}} = \sum_{|\mathbf{k}|=1}^m \gamma_{\mathbf{k}} C^{\mathbf{k}}_{\mathbf{j}} \;, \quad \text{for each } \mathbf{j} \; \text{where } |\mathbf{j}|\leq m \,,
\end{align*}
```
were we have introduced multi-index scalars
```math
\begin{align*}
    C^{\mathbf{u}}_{\mathbf{v}} &= C^{u_1}_{v_1}C^{u_2}_{v_2} \dotsm C^{u_N}_{v_N},
\end{align*}
```
with each element defined as
```math
\begin{align*}
C^{a}_{b} &=
\begin{cases}
  \frac{a!}{(a-b)!b!}, & a \geq b \\
  0, & a \lt b
\end{cases}\;.
\end{align*}
```

Using the specific construction of $\bar{\mathbf{φ}}$ described above, it can be shown [15] that for every *deterministic* initial condition, i.e., $\mathbf{n}(t_0) = \mathbf{n}_0$ with probability one, we will have
```math
\begin{align*}
    \mathbf{μ}(t_0) = \mathbf{ν}(t_0) &\implies \left. \frac{d \mathbf{μ}}{dt} \right\rvert_{t=t_0} = \left. \frac{d \mathbf{ν}}{dt} \right\rvert_{t=t_0} \\
    &\implies \left. \frac{d^2\mathbf{μ}}{dt^2} \right\rvert_{t=t_0} = \left. \frac{d^2 \mathbf{ν}}{dt^2} \right\rvert_{t=t_0} + \mathbf{ϵ}(\mathbf{n}_0) \;,
\end{align*}
```
where all elements of $\mathbf{ϵ}(\mathbf{n}_0)$ are zero except the ones corresponding to $m^{\text{th}}$ order raw moments—these elements are second order polynomials in $\mathbf{n}_0$. Note, however, that these results hold only for mass-action systems containing no higher than second order chemical reactions. While the derivative matching MA can be applied in the same way to systems containing higher order polynomial and non-polynomial propensity functions, it has not been rigorously analysed in such scenarios, where, naturally, we expect significantly larger approximation errors.

## Conditional closures

Discussion of [16].

### Conditional normal closure

### Conditional derivative matching

## References

[1]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Approximation and inference methods for stochastic biochemical kinetics - a tutorial review", Journal of Physics A: Mathematical and Theoretical 50, 093001 (2017). https://doi.org/10.1088/1751-8121/aa54d9

[2]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Comparison of different moment-closure approximations for stochastic chemical kinetics", The Journal of Chemical Physics 143, 185101 (2015). https://doi.org/10.1063/1.4934990

[3]: E. Lakatos, A. Ale, P. D. W. Kirk, and M. P. H. Stumpf, "Multivariate moment closure techniques for stochastic kinetic models", The Journal of Chemical Physics 143, 094107 (2015). https://doi.org/10.1063/1.4929837

[4]: J.  Hespanha,  "Moment  closure  for  biochemical  networks",  in  2008  3rd  International Symposium on Communications, Control and Signal Processing (Mar. 2008), pp. 142–147. https://doi.org/10.1109/ISCCSP.2008.4537208

[5]: L. A. Goodman, "Population Growth of the Sexes", Biometrics9, Publisher: [Wiley, International Biometric Society], 212–225 (1953). https://doi.org/10.2307/3001852

[6]: P. Whittle, "On the use of the normal approximation in the treatment of stochastic processes", Journal of the Royal Statistical Society: Series B (Methodological) 19, 268–281 (1957). https://doi.org/10.1111/j.2517-6161.1957.tb00263.x

[7]: N. Balakrishnan, N. L. Johnson, and S. Kotz, “A note on relationships between moments, central moments and cumulants from multivariate distributions”, Statistics & Probability Letters 39, 49–54 (1998). https://doi.org/10.1016/S0167-7152(98)00027-3

[8]: I. Nasell, "An extension of the moment closure method", Theoretical Population Biology 64, 233–239 (2003). https://doi.org/10.1016/S0040-5809(03)00074-1

[9]: M. J. Keeling, "Multiplicative Moments and Measures of Persistence in Ecology", Journal of Theoretical Biology 205, 269–281 (2000). https://doi.org/10.1006/jtbi.2000.2066

[10]: E. L. Crow and K. Shimizu, eds., Lognormal Distributions: Theory and Applications (Marcel Dekker, 1988).

[11]: N. L. Johnson, S. Kotz, and N. Balakrishnan, Discrete Multivariate Distributions (Wiley, Feb. 1997).

[12]: A. M. Mathal and P. G. Moschopoulos, "A form of multivariate gamma distribution", Annals of the Institute of Statistical Mathematics 44, 97–106 (1992). https://doi.org/10.1007/BF00048672

[13]: E. Furman, "On a multivariate gamma distribution", Statistics & Probability Letters 78, 2353–2360 (2008). https://doi.org/10.1016/j.spl.2008.02.012

[14]: A.  Singh  and  J.  P.  Hespanha,  "Lognormal  Moment  Closures  for  Biochemical  Reactions", in Proceedings of the 45th IEEE Conference on Decision and Control, ISSN:0191-2216 (Dec. 2006), pp. 2063–2068.

[15]: A. Singh and J. P. Hespanha, "Approximate Moment Dynamics for Chemically Reacting Systems", IEEE Transactions on Automatic Control 56, 414–418 (2011).

[16]: M. Soltani, C. A. Vargas-Garcia, and A. Singh, "Conditional Moment Closure Schemes for Studying Stochastic Dynamics of Genetic Circuits", IEEE Transactions on Biomedical Circuits and Systems 9, 518–526 (2015).
