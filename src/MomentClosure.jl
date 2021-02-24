module MomentClosure

import Catalyst: species, params, reactions, speciesmap, paramsmap, numspecies, numreactions, numparams,
                 substoichmat, prodstoichmat

using ModelingToolkit
using ModelingToolkit: value
using SciMLBase: EnsembleSolution

using SymbolicUtils
using SymbolicUtils.Rewriters: Chain, RestartedChain, PassThrough, Prewalk, Postwalk, Fixpoint
using SymbolicUtils: @rule, @acrule, @ordered_acrule, isnotflat, flatten_term, istree,
                     needs_sorting, sort_args, is_literal_number, hasrepeats, merge_repeats,
                     _iszero, pow, one, _isone, zero, symtype, operation, arguments

using Combinatorics: permutations
using TupleTools: sort
using Cumulants

export cumulants_to_raw_moments, cumulants_to_central_moments,
       raw_to_central_moments, central_to_raw_moments,
       generate_central_moment_eqs, generate_raw_moment_eqs, bernoulli_moment_eqs,
       get_S_mat, propensities, ReactionSystemMod,
       species, params, speciesmap, paramsmap, numspecies, numreactions, numparams,
       @parameters, @variables, # ModelingToolkit variables needed for model initialisation
       clean_expr, expand_mod, simplify, # Symbolic manipulation tools (simplify borrowed from ModelingToolkit & SymbolicUtils)
       moment_closure, deterministic_IC,
	   format_moment_eqs, format_closure,
	   sample_raw_moments, sample_central_moments, sample_cumulants

include("symbolic.jl")
include("moment_convert.jl")
include("reaction_systems.jl")
include("stochastic_stoichiometry.jl")
include("central_moment_equations.jl")
include("raw_moment_equations.jl")
include("utils.jl")
include("bernoulli.jl")
include("closure_methods/closure.jl")
include("closure_methods/zero_closure.jl")
include("closure_methods/normal_closure.jl")
include("closure_methods/log_normal_closure.jl")
include("closure_methods/poisson_closure.jl")
include("closure_methods/gamma_closure.jl")
include("closure_methods/derivative_matching.jl")
include("closure_methods/conditional_gaussian.jl")
include("closure_methods/conditional_derivative_matching.jl")

end