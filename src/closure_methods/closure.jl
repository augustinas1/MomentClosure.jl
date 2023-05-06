function close_eqs(sys::MomentEquations, closure_exp::OrderedDict,
                   closure::OrderedDict, polynorm::Bool)

    closed_eqs = Equation[]
    for eq in sys.odes.eqs

        closed_rhs = substitute(eq.rhs, closure_exp)
        # apply binomial expansion on the expressions
        if polynorm # depending on the functional form
            closed_rhs = expand(closed_rhs)
        else
            closed_rhs = expand_mod(closed_rhs)
        end
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end

    iv = get_iv(sys.odes)
    ps = get_ps(sys.odes)
    
    vars = states(sys.odes)[1:(length(sys.iter_1)+length(sys.iter_m))]

    odename = Symbol(nameof(sys), "_closed")
    odes = ODESystem(closed_eqs, iv, vars, ps; name=odename)

    ClosedMomentEquations(odes, closure, sys)

end

"""
    moment_closure(sys::MomentEquations, closure::String, binary_vars::Array{Int,1}=Int[])

Given `MomentEquations`, apply the specified moment closure approximation and return
the [`ClosedMomentEquations`](@ref).

The supported `closure` options are:
- [`"zero"`](@ref zero_closure)
- [`"normal"`](@ref normal_closure)
- [`"log-normal"`](@ref log-normal_closure)
- [`"poisson"`](@ref poisson_closure)
- [`"gamma"`](@ref gamma_closure)
- [`"derivative matching"`](@ref derivative_matching)
- [`"conditional gaussian"`](@ref conditional_gaussian_closure)
- [`"conditional derivative matching"`](@ref conditional_derivative_matching)

# Notes
- `binary_vars` *must* be specified for conditional closures as an array of indices of all species
  (as in [`Catalyst.speciesmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.speciesmap)) 
  which molecule number is a Bernoulli variable. Although not necessary
  for other closures, specifying `binary_vars` is recommended as the properties of Bernoulli variables
  will be used to remove the redundant moment equations and simplify the expressions, which can
  significantly improve numerical stability.
"""
function moment_closure(sys::MomentEquations, closure::String, binary_vars::Array{Int,1}=Int[])

    if closure == "zero"
        zero_closure(sys, binary_vars)
    elseif closure == "normal"
        normal_closure(sys, binary_vars)
    elseif closure == "log-normal"
        log_normal_closure(sys, binary_vars)
    elseif closure == "poisson"
        poisson_closure(sys, binary_vars)
    elseif closure == "gamma"
        gamma_closure(sys, binary_vars)
    elseif closure == "derivative matching"
        derivative_matching(sys, binary_vars)
    elseif closure == "conditional gaussian"
        conditional_gaussian_closure(sys, binary_vars)
    elseif closure == "conditional derivative matching"
        conditional_derivative_matching(sys, binary_vars)
    else
        error(closure*" closure does not exist")
    end

end
