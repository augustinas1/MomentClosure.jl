module MomentClosure

using Catalyst
import Catalyst: species, params, reactions, speciesmap, paramsmap, numspecies,
				 numreactions, numparams, substoichmat, prodstoichmat, netstoichmat,
				 ReactionSystem, get_eqs, get_states, get_iv, get_ps

using ModelingToolkit
using Symbolics: value, var_from_nested_derivative, map_subscripts, hessian, gradient, setmetadata, scalarize
import Symbolics: degree

using SciMLBase, SciMLBase.EnsembleAnalysis
using DiffEqJump
using Random
using Distributions: Geometric

# TODO: remove unnecessary imports for SymbolicUtils
using SymbolicUtils.Rewriters: Chain, PassThrough, Prewalk, Fixpoint
using SymbolicUtils: Symbolic, Term, Real, expand, simplify, operation,
					 arguments, @rule, @acrule, isnotflat, flatten_term,
					 istree, FnType

using DataStructures: OrderedSet, OrderedDict
using TupleTools: sort
using Combinatorics
using Cumulants
using Latexify

using DocStringExtensions

export generate_central_moment_eqs, generate_raw_moment_eqs, bernoulli_moment_eqs,
       propensities, ReactionSystemMod, species, params, speciesmap, paramsmap,
	   numspecies, numreactions, numparams, netstoichmat,
       moment_closure, deterministic_IC, JumpProblem,
	   get_raw_moments, get_central_moments, get_cumulants, get_moments_FSP,
	   linear_mapping_approximation

# needed for ReactionSystemMod definition
export @parameters, @variables

include("reaction_systems.jl")
include("moment_equations.jl")
include("symbolic.jl")
include("moment_convert.jl")
include("stochastic_stoichiometry.jl")
include("central_moment_equations.jl")
include("central_moment_equations_SDE.jl")
include("raw_moment_equations.jl")
include("raw_moment_equations_SDE.jl")
include("raw_moment_equations_JumpDiff.jl")
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
include("closure_methods/linear_mapping_approximation.jl")
include("stochastic_simulation.jl")

end
