"""
$(TYPEDEF)

Closed moment equations and the corresponding closure functions.

# Fields
$(FIELDS)
"""
struct ClosedMomentEquations <: MomentEquations
    """[`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/)
    consisting of the time-evolution equations of *closed* moments."""
    odes::ODESystem
    """Dictionary of moment closure functions for each higher order moment."""
    closure::Dict
    """Original raw or central moment equations (before closure was applied)."""
    open_eqs::MomentEquations
end

function close_eqs(sys::MomentEquations, closure::Dict, closure_symbolic::Dict)

    closed_eqs = []
    for eq in sys.odes.eqs
        closed_rhs = substitute(eq.rhs, closure)
        closed_rhs = simplify(closed_rhs)
        push!(closed_eqs, Equation(eq.lhs, closed_rhs))
    end

    ClosedMomentEquations(ODESystem(closed_eqs), closure_symbolic, sys)

end


function moment_closure(sys::MomentEquations, closure_name::String,
    binary_vars::Vector=[]; clean=true)

    if closure_name == "zero"
        zero_closure(sys)
    elseif closure_name == "normal"
        normal_closure(sys)
    elseif closure_name == "log-normal"
        log_normal_closure(sys)
    elseif closure_name == "poisson"
        poisson_closure(sys)
    elseif closure_name == "gamma"
        gamma_closure(sys, clean)
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
