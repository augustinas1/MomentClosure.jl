# MomentClosure.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://augustinas1.github.io/MomentClosure.jl/dev/)
[![Build Status](https://github.com/augustinas1/MomentClosure.jl/workflows/CI/badge.svg)](https://github.com/augustinas1/MomentClosure.jl/actions)
[![Coverage](https://codecov.io/gh/augustinas1/MomentClosure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/augustinas1/MomentClosure.jl)

MomentClosure.jl is a tool to automatically obtain time-evolution equations of moments up to an arbitrary order for virtually any chemical reaction network implementing a wide array of moment closure approximations commonly used in stochastic biochemical kinetics [[1]](#1). MomentClosure is (attempted to be) fairly well-integrated within the broader Julia ecosystem utilising a number of familiar packages:
- MomentClosure can be immediately applied to reaction network models defined using [Catalyst](https://github.com/SciML/Catalyst.jl/issues/22).
- Moment equations are generated as a [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) [`ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) (with some extra help from [SymbolicUtils](https://github.com/JuliaSymbolics/SymbolicUtils.jl)).
- The resulting `ODESystem` can be solved using any [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl/) ODE solvers, enabling further study of the system using [parameter estimation](https://diffeq.sciml.ai/stable/analysis/parameter_estimation/), [sensitivity analysis](https://diffeq.sciml.ai/stable/analysis/sensitivity/) and [bifurcation analysis](https://diffeq.sciml.ai/stable/analysis/bifurcation/) tools (have no example/test for these yet!).

## Tutorials and documentation

Please [see the documentation](https://augustinas1.github.io/MomentClosure.jl/dev/) for information on using the package, theory behind it and in-depth examples (yet to be added).

## Features
- Chemical reaction networks containing any number of molecular species and reactions with any type of *smooth* propensity functions are supported.
- Models can be defined using [Catalyst](https://github.com/SciML/Catalyst.jl/issues/22) as [`ModelingToolkit.ReactionSystem`](https://catalyst.sciml.ai/dev/api/catalyst_api/#ModelingToolkit.ReactionSystem). Alternatively, built-in functionality (heavily based on Catalyst) can be used for model initialisation.
- Added (experimental and limited) support for reaction networks involving reaction products that are geometrically distributed random variables. An example of such network is an autoregulatory gene network with bursty protein production where the burst size follows a geometric distribution. To our knowledge, such models cannot be defined using Catalyst and hence modelled using the broader [SciML framework](https://github.com/SciML/). For this reason, we also include an implementation of the Gillespie's Stochastic Simulation Algorithm [[2]](#2) generally applicable to such reaction networks (yet to be added).
- Equations describing the time evolution of means and central moments of the number of molecules of each species in the system can be generated up to arbitrary order [[3, 4]](#3). Note that non-polynomial propensity functions are Taylor expanded to a specified order. Raw moment equations can also be generated for mass-action systems (where all propensity functions are polynomials).
- Currently implemented moment closure approximations include:
	- zero (central-moment neglect) closure [[3, 5]](#3)
	- normal closure [[5]](#5)
	- poisson closure [[5]](#5)
	- log-normal closure [[6]](#6)
	- gamma closure [[6]](#6)
	- derivative matching [[7]](#7)
	- conditional gaussian closure [[8]](#8)
	- conditional derivative matching [[8]](#8)
- Moment equations are constructed as a [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/) that can be solved using any [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl/) ODE solver. Moreover, [parameter estimation](https://diffeq.sciml.ai/stable/analysis/parameter_estimation/), [sensitivity analysis](https://diffeq.sciml.ai/stable/analysis/sensitivity/) and [bifurcation analysis](https://diffeq.sciml.ai/stable/analysis/bifurcation/) tools can be applied to further study the resulting system of equations.
- [Latexify](https://github.com/korsbo/Latexify.jl) can be used to generate LaTeX expressions of the corresponding moment equations.

## References

<a id="1">[1]</a> D. Schnoerr, G. Sanguinetti, and R. Grima, "Approximation and inference methods for stochastic biochemical kinetics - a tutorial review", Journal of Physics A: Mathematical and Theoretical 50, 093001 (2017). https://doi.org/10.1088/1751-8121/aa54d9

<a id="2">[2]</a>: D. T. Gillespie, "A general method for numerically simulating the stochastic time evolution of coupled chemical reactions", Journal of Computational Physics 22, 403-434 (1976). https://doi.org/10.1016/0021-9991(76)90041-3

<a id="3">[3]</a>: A. Ale, P. Kirk, and M. P. H. Stumpf, "A general moment expansion method for stochastic kinetic models", The Journal of Chemical Physics 138, 174101 (2013). https://doi.org/10.1063/1.4802475

<a id="4">[4]</a>: C. H. Lee, "A Moment Closure Method for Stochastic Chemical Reaction Networks with General Kinetics", MATCH Communications in Mathematical and in Computer Chemistry 70, 785-800 (2013). https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_785-800.pdf

<a id="5">[5]</a>: D. Schnoerr, G. Sanguinetti, and R. Grima, "Comparison of different moment-closure approximations for stochastic chemical kinetics", The Journal of Chemical Physics 143, 185101 (2015). https://doi.org/10.1063/1.4934990

<a id="6">[6]</a>: E. Lakatos, A. Ale, P. D. W. Kirk, and M. P. H. Stumpf, "Multivariate moment closure techniques for stochastic kinetic models", The Journal of Chemical Physics 143, 094107 (2015). https://doi.org/10.1063/1.4929837

<a id="7">[7]</a>: A. Singh and J. P. Hespanha, "Lognormal Moment Closures for Biochemical Reactions", in Proceedings of the 45th IEEE Conference on Decision and Control, ISSN: 0191-2216 (Dec. 2006), pp. 2063-2068. https://doi.org/10.1109/CDC.2006.376994

<a id="8">[8]</a>: M. Soltani, C. A. Vargas-Garcia, and A. Singh, "Conditional Moment Closure Schemes for Studying Stochastic Dynamics of Genetic Circuits", IEEE Transactions on Biomedical Circuits and Systems 9, 518-526 (2015). https://doi.org/10.1109/TBCAS.2015.2453158
