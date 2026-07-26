module MomentClosure

import ModelingToolkitBase
using ModelingToolkitBase: ODESystem, System, complete, equations, get_eqs, get_iv,
    get_noise_eqs, get_ps, parameters, unknowns
using SymbolicIndexingInterface: getname
import Catalyst
using Catalyst: ReactionSystem, default_t, jumpratelaw, numreactions, numspecies, reactions,
    species, speciesmap, substoichmat
import SciMLBase
using SciMLBase: EnsembleSolution, NullParameters, ODEProblem, ODESolution
using SciMLBase.EnsembleAnalysis: componentwise_vectors_timestep
using Distributions: Geometric

import Symbolics
using Symbolics: Differential, Equation, Num, @variables, expand_derivatives, factors,
    get_variables, scalarize, setmetadata, value
import SymbolicUtils.Rewriters: Chain, Fixpoint, PassThrough, Prewalk
using SymbolicUtils: @acrule, @rule, BasicSymbolic, arguments, expand, isadd, isdiv, ismul,
    isterm, operation, simplify, substitute
using TermInterface: iscall

using DataStructures: OrderedDict
using Combinatorics: multiset_permutations, partitions
using Cumulants: cumulants, naivecumulant, naivemoment
import Latexify
using Latexify: @latexrecipe
using StatsBase: moment

using DocStringExtensions: FIELDS, TYPEDEF

export generate_central_moment_eqs, generate_raw_moment_eqs, bernoulli_moment_eqs,
    propensities, get_stoichiometry, moment_closure, deterministic_IC,
    get_raw_moments, get_central_moments, get_cumulants, get_moments_FSP,
    linear_mapping_approximation, get_odes, get_closure

include("reaction_systems.jl")
include("moment_equations.jl")
include("symbolic.jl")
include("moment_convert.jl")
include("stochastic_stoichiometry.jl")
include("central_moment_equations.jl")
include("central_moment_equations_SDE.jl")
include("raw_moment_equations.jl")
include("raw_moment_equations_SDE.jl")
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

end
