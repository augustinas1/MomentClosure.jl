function close_eqs(sys::Union{RawMomentEquations, CentralMomentEquations},
                   closure, closure_symbolic)

    closed_eqs = []
    for eq in sys.odes.eqs
        closed_rhs = substitute(eq.rhs, closure)
        closed_rhs = simplify(closed_rhs)
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end
    ODESystem(closed_eqs), closure, closure_symbolic

end

function moment_closure(sys::Union{RawMomentEquations, CentralMomentEquations}, closure_name::String)

    if closure_name == "zero"
        zero_closure(sys)
    else
        error("this closure does not exist")
    end

end
