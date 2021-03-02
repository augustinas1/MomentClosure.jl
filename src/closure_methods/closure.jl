function close_eqs(sys::MomentEquations,
                   closure, closure_symbolic)

    closed_eqs = []
    for eq in sys.odes.eqs
        closed_rhs = substitute(eq.rhs, closure)
        closed_rhs = simplify(closed_rhs)
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end
    ODESystem(closed_eqs), closure_symbolic

end

function moment_closure(sys::MomentEquations,
        closure_name::String, binary_vars::Vector=[])

    if closure_name == "zero"
        zero_closure(sys)
    elseif closure_name == "normal"
        normal_closure(sys)
    elseif closure_name == "log-normal"
        log_normal_closure(sys)
    elseif closure_name == "poisson"
        poisson_closure(sys)
    elseif closure_name == "gamma"
        gamma_closure(sys)
    elseif closure_name == "derivative matching"
        derivative_matching(sys)
    elseif closure_name == "conditional gaussian"
        conditional_gaussian_closure(sys, binary_vars)
    elseif closure_name == "conditional derivative matching"
        conditional_derivative_matching(sys, binary_vars)
    else
        error(closure_name*" closure does not exist")
    end

end
