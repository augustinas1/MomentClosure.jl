# MomentClosure

MomentClosure.jl is a tool to automatically obtain time-evolution equations of moments up to an arbitrary order for virtually any chemical reaction network or system of stochastic differential equations (SDEs), implementing a wide array of moment closure approximations commonly used in stochastic biochemical kinetics [1]. MomentClosure is (attempted to be) fairly well-integrated within the broader Julia ecosystem utilising a number of familiar packages:
- MomentClosure can be immediately applied to reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl) and SDE systems built with [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl).
- Moment equations are generated as a [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) [`ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) (with some extra help from [Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) and [SymbolicUtils](https://github.com/JuliaSymbolics/SymbolicUtils.jl)).
- The resulting `ODESystem` can be solved using any [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl/) ODE solvers, enabling further study of the system using [parameter estimation](https://diffeq.sciml.ai/stable/analysis/parameter_estimation/), [sensitivity analysis](https://diffeq.sciml.ai/stable/analysis/sensitivity/) and [bifurcation analysis](https://diffeq.sciml.ai/stable/analysis/bifurcation/) tools.

## Features
* Chemical reaction networks containing any number of molecular species and reactions with any type of *smooth* propensity functions are supported. Models can be defined using [Catalyst](https://github.com/SciML/Catalyst.jl/issues/22) as [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/dev/api/catalyst_api/#ModelingToolkit.ReactionSystem). Alternatively, built-in functionality (heavily based on Catalyst) can be used for model initialisation.
* Added (limited) support for reaction networks involving reaction products that are geometrically distributed random variables. An example of such network is an autoregulatory gene network with bursty protein production where the burst size follows a geometric distribution, see [this example](@ref geometric-and-conditional) for more details. Such models *currently* cannot be defined using Catalyst and hence modelled using the broader [SciML framework](https://github.com/SciML/).
* Equations describing the time evolution of means and central moments of the number of molecules of each species in the system can be generated up to arbitrary order [2, 3]. Note that non-polynomial propensity functions are Taylor expanded to a specified order. Raw moment equations can also be generated for mass-action systems (where all propensity functions are polynomials). 
* SDE systems defined as [`ModelingToolkit.SDESystem`](https://mtk.sciml.ai/stable/systems/SDESystem/#ModelingToolkit.SDESystem) are supported (big thanks to @FHoltorf). Similarly to reaction networks with non-polynomial propensities, non-polynomial drift and diffusion coefficients are Taylor expanded.
* Currently implemented moment closure approximations include:
  - zero (central-moment neglect) closure [4, 5]
  - normal closure [4]
  - poisson closure [4]
  - log-normal closure [5]
  - gamma closure [5]
  - derivative matching [6]
  - conditional gaussian closure [7]
  - conditional derivative matching [7]
  - linear mapping approximation [8]
* Moment equations are constructed as a [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) that can be solved using any [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl/) ODE solver. Moreover, [parameter estimation](https://diffeq.sciml.ai/stable/analysis/parameter_estimation/), [sensitivity analysis](https://diffeq.sciml.ai/stable/analysis/sensitivity/) and [bifurcation analysis](https://diffeq.sciml.ai/stable/analysis/bifurcation/) tools can be applied to further study the resulting system of equations.
- [Latexify](https://github.com/korsbo/Latexify.jl) can be used to generate LaTeX expressions of the corresponding moment equations.


## Installation
MomentClosure can be installed through the Julia package manager:

```julia
]add https://github.com/augustinas1/MomentClosure.jl
using MomentClosure
```

Note the use of full link to github repo as the package is currently unregistered.

## Citation

If you use MomentClosure in your work, please cite [our paper](https://arxiv.org/abs/2105.05475):
```
@article{MomentClosure2021,
    author = {Sukys, Augustinas and Grima, Ramon},
    title = "{MomentClosure.jl: automated moment closure approximations in Julia}",
    journal = {Bioinformatics},
    volume = {38},
    number = {1},
    pages = {289-290},
    year = {2021},
    month = {06},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btab469},
    url = {https://doi.org/10.1093/bioinformatics/btab469},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/38/1/289/41891091/btab469.pdf},
}
```

## References

[1]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Approximation and inference methods for stochastic biochemical kinetics - a tutorial review", Journal of Physics A: Mathematical and Theoretical 50, 093001 (2017). [https://doi.org/10.1088/1751-8121/aa54d9](https://doi.org/10.1088/1751-8121/aa54d9)

[2]: A. Ale, P. Kirk, and M. P. H. Stumpf, "A general moment expansion method for stochastic kinetic models", The Journal of Chemical Physics 138, 174101 (2013). [https://doi.org/10.1063/1.4802475](https://doi.org/10.1063/1.4802475)

[3]: C. H. Lee, "A Moment Closure Method for Stochastic Chemical Reaction Networks with General Kinetics", MATCH Communications in Mathematical and in Computer Chemistry 70, 785-800 (2013). [https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_785-800.pdf](https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_785-800.pdf)

[4]: D. Schnoerr, G. Sanguinetti, and R. Grima, "Comparison of different moment-closure approximations for stochastic chemical kinetics", The Journal of Chemical Physics 143, 185101 (2015). [https://doi.org/10.1063/1.4934990](https://doi.org/10.1063/1.4934990)

[5]: E. Lakatos, A. Ale, P. D. W. Kirk, and M. P. H. Stumpf, "Multivariate moment closure techniques for stochastic kinetic models", The Journal of Chemical Physics 143, 094107 (2015). [https://doi.org/10.1063/1.4929837](https://doi.org/10.1063/1.4929837)

[6]: A. Singh and J. P. Hespanha, "Lognormal Moment Closures for Biochemical Reactions", in Proceedings of the 45th IEEE Conference on Decision and Control, ISSN: 0191-2216 (Dec. 2006), pp. 2063-2068. [https://doi.org/10.1109/CDC.2006.376994](https://doi.org/10.1109/CDC.2006.376994)

[7]: M. Soltani, C. A. Vargas-Garcia, and A. Singh, "Conditional Moment Closure Schemes for Studying Stochastic Dynamics of Genetic Circuits", IEEE Transactions on Biomedical Circuits and Systems 9, 518-526 (2015). [https://doi.org/10.1109/TBCAS.2015.2453158](https://doi.org/10.1109/TBCAS.2015.2453158)

[8]: Z. Cao and R. Grima, "Linear mapping approximation of gene regulatory networks with stochastic dynamics", Nature Communications 9, 3305 (2018). [https://doi.org/10.1038/s41467-018-05822-0](https://doi.org/10.1038/s41467-018-05822-0)
