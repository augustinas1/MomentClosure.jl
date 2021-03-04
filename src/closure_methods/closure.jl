function close_eqs(sys::MomentEquations, closure::Dict, closure_symbolic::Dict)

    closed_eqs = Equation[]
    for eq in sys.odes.eqs
        closed_rhs = substitute(eq.rhs, closure)
        closed_rhs = simplify(closed_rhs)
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end

    iv = sys.odes.iv
    ps = sys.odes.ps

    vars = extract_variables(closed_eqs, sys.N, sys.q_order)
    #vars = extract_variables(closed_eqs, ps)
    odes = ODESystem(closed_eqs, iv, vars, ps)

    ClosedMomentEquations(odes, closure_symbolic, sys)

end

"""
    moment_closure(sys::MomentEquations, closure::String, binary_vars::Vector=[]; clean=true)

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
- `binary_vars` is required in case of conditional closures. It must be
  specified by denoting in an array the indices of all species whose molecule
  number is a Bernoulli (binary) variable.
- `clean` argument is only relevant for gamma closure. Namely, `clean=true`
  implies that the closure functions and the resulting closed moment equations
  will be heavily simplified in order to preserve numerical stability. However,
  this can extremely computationally expensive due to complicated functional
  form of the higher order moments under gamma closure. Setting `clean=false`
  is generally much faster as the symbolic expressions are not simplified
  internally, but it can lead to severe numerical instabilities.
"""
function moment_closure(sys::MomentEquations, closure::String,
    binary_vars::Vector=[]; clean=true)

    if closure == "zero"
        zero_closure(sys)
    elseif closure == "normal"
        normal_closure(sys)
    elseif closure == "log-normal"
        log_normal_closure(sys)
    elseif closure == "poisson"
        poisson_closure(sys)
    elseif closure == "gamma"
        gamma_closure(sys, clean)
    elseif closure == "derivative matching"
        derivative_matching(sys)
    elseif closure == "conditional gaussian"
        conditional_gaussian_closure(sys, binary_vars)
    elseif closure == "conditional derivative matching"
        conditional_derivative_matching(sys, binary_vars)
    else
        error(closure*" closure does not exist")
    end

end
