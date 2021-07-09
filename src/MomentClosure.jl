module MomentClosure

import Catalyst: species, params, reactions, speciesmap, paramsmap, numspecies,
				 numreactions, numparams, substoichmat, prodstoichmat

using ModelingToolkit
using Symbolics: value, var_from_nested_derivative, map_subscripts

using SciMLBase, SciMLBase.EnsembleAnalysis
using DiffEqJump
using Random
using Distributions: Geometric

using SymbolicUtils.Rewriters: Chain, PassThrough, Prewalk, Fixpoint
using SymbolicUtils: Symbolic, Term, Real, expand, simplify, operation,
					 arguments, @rule, @acrule, isnotflat, flatten_term,
					 istree, FnType

using OrderedCollections: OrderedDict
using Combinatorics
using TupleTools: sort
using Cumulants
using Latexify

using DocStringExtensions

export generate_central_moment_eqs, generate_raw_moment_eqs, bernoulli_moment_eqs,
       get_S_mat, propensities, ReactionSystemMod,
       species, params, speciesmap, paramsmap, numspecies, numreactions, numparams,
       moment_closure, deterministic_IC, JumpProblem,
	   get_raw_moments, get_central_moments, get_cumulants

# reexporting from ModelingToolkit & Symbolics
# needed for ReactionSystemMod definition
export @parameters, @variables

include("reaction_systems.jl")
include("moment_equations.jl")
include("symbolic.jl")
include("moment_convert.jl")
include("stochastic_stoichiometry.jl")
include("central_moment_equations.jl")
include("raw_moment_equations.jl")
include("closure_methods/closure.jl")
include("bernoulli.jl")
include("utils.jl")
include("closure_methods/zero_closure.jl")
include("closure_methods/normal_closure.jl")
include("closure_methods/log_normal_closure.jl")
include("closure_methods/poisson_closure.jl")
include("closure_methods/gamma_closure.jl")
include("closure_methods/derivative_matching.jl")
include("closure_methods/conditional_gaussian.jl")
include("closure_methods/conditional_derivative_matching.jl")
include("stochastic_simulation.jl")


end
