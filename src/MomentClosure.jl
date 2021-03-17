module MomentClosure

import Catalyst: species, params, reactions, speciesmap, paramsmap, numspecies,
				 numreactions, numparams, substoichmat, prodstoichmat

using ModelingToolkit
using Symbolics: value, var_from_nested_derivative

using SciMLBase, SciMLBase.EnsembleAnalysis

using SymbolicUtils.Rewriters: Chain, PassThrough, Prewalk, Fixpoint
using SymbolicUtils: polynormalize, simplify, operation, arguments, istree
using SymbolicUtils: @rule, @acrule, isnotflat, flatten_term

using OrderedCollections: OrderedDict
using Combinatorics: permutations
using TupleTools: sort
using Cumulants
using Latexify

using DocStringExtensions

export generate_central_moment_eqs, generate_raw_moment_eqs, bernoulli_moment_eqs,
       get_S_mat, propensities, ReactionSystemMod,
       species, params, speciesmap, paramsmap, numspecies, numreactions, numparams,
       moment_closure, deterministic_IC, ODEProblem,
	   sample_raw_moments, sample_central_moments, sample_cumulants

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



end
