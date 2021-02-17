function close_eqs(sys::Union{RawMomentEquations, CentralMomentEquations},
                   closure, closure_symbolic)

    closed_eqs = []
    for eq in sys.odes.eqs
        closed_rhs = substitute(eq.rhs, closure)
        closed_rhs = simplify(closed_rhs)
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end
    ODESystem(closed_eqs), closure_symbolic

end

function moment_closure(sys::Union{RawMomentEquations, CentralMomentEquations}, closure_name::String)

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
    elseif closure_name == "derivative_matching"
        derivative_matching(sys)
    else
        error(closure_name*" closure does not exist")
    end

end
