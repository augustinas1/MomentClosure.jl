function bernoulli_iter_redundant(sys, binary_vars::Array{Int,1})

    redundant_iter_sub = Dict()

    # starting from N+2 index as 0th and 1st order terms stay the same
    for iter in sys.iter_all[sys.N+2:end]

        iter_temp = [iter...] # need to transform into array as Tuple is immutable
        for ind in binary_vars
            if iter[ind] > 1
                iter_temp[ind] = 1
            end
        end
        iter_temp = Tuple(iter_temp)

        if iter != iter_temp
            redundant_iter_sub[iter] = iter_temp
        end

    end

    redundant_iter_sub

end


function bernoulli_reduce(sys::MomentEquations, binary_vars::Array{Int,1})

    N = sys.N

    redundant_iter_sub = bernoulli_iter_redundant(sys, binary_vars)
    redundant_iter = keys(redundant_iter_sub)

    redundant_eqs = Int[]
    for ind in binary_vars
        for (i, iter) in enumerate(sys.iter_m)
            if iter[ind] > 1
                push!(redundant_eqs, i+N)
            end
        end
    end

    if sys isa RawMomentEquations
        μ = copy(sys.μ)
        iter_sub = [μ[key] => μ[val] for (key, val) in redundant_iter_sub]
    else
        μ = define_μ(N, sys.q_order)
        μ_redundant_sub = [μ[key] => μ[val] for (key, val) in redundant_iter_sub]

        clean_iter = setdiff(sys.iter_all[N+2:end], redundant_iter)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        μ_clean_sub = Dict([μ[iter] => central_to_raw[iter] for iter in clean_iter])

        raw_to_central = raw_to_central_moments(N, sys.q_order)
        iter_sub = Dict()
        for iter in redundant_iter
            M_temp = raw_to_central[iter]
            M_temp = simplify(substitute(M_temp, μ_redundant_sub))
            M_temp = substitute(M_temp, μ_clean_sub)
            M_temp = simplify(M_temp)
            iter_sub[sys.M[iter]] = M_temp
        end

    end

    redundant_iter, redundant_eqs, iter_sub

end

"""
    bernoulli_moment_eqs(sys::MomentEquations, binary_vars::Array{Int,1})

Given `MomentEquations` and an array of indices specifying the species which molecule numbers are
binary variables (either 0 or 1), apply identities of Bernoulli variables to remove the redundant
ODEs and return the *cleaned up* `MomentEquations`. See [here](@ref geometric-and-conditional) for
example usage.
"""
function bernoulli_moment_eqs(sys::MomentEquations, binary_vars::Array{Int,1})

    # Construct moment equations removing the redundant ones
    # noting the properties of the Bernoulli variables in the system

    redundant_iter, redundant_eqs, iter_sub = bernoulli_reduce(sys, binary_vars)

    # construct the cleaned moment equations
    clean_eqs = Equation[]
    for (i, eq) in enumerate(get_eqs(sys))
        if !(i in redundant_eqs)
            clean_rhs = substitute(eq.rhs, iter_sub)
            clean_rhs = expand(clean_rhs)
            push!(clean_eqs, Equation(eq.lhs, clean_rhs))
        end
    end

    # construct a new Raw/Central/MomentEquations system saving only the clean iterators
    clean_iter = setdiff(sys.iter_all, redundant_iter)
    iter_m = filter(x -> 2 <= sum(x) <= sys.m_order, clean_iter)
    iter_q = filter(x -> sys.m_order < sum(x) <= sys.q_order, clean_iter)

    field_values = [getfield(sys, field) for field in fieldnames(typeof(sys))]
    ind_iter_m = findall(x -> x==:iter_m, fieldnames(typeof(sys)))[1]
    ind_iter_q = findall(x -> x==:iter_q, fieldnames(typeof(sys)))[1]
    ind_iter_all = findall(x -> x==:iter_all, fieldnames(typeof(sys)))[1]

    field_values[ind_iter_m] = iter_m       #sys.iter_m
    field_values[ind_iter_q] = iter_q       #sys.iter_q
    field_values[ind_iter_all] = clean_iter #sys.iter_all

    ## fixing ODE system to preserve consistent ordering of parameters
    iv = get_iv(get_odes(sys))
    ps = get_ps(get_odes(sys))
    
    if sys isa RawMomentEquations
        vars = extract_variables(clean_eqs, sys.μ)
    else
        vars = extract_variables(clean_eqs, sys.μ, sys.M)
    end
   
    odename = Symbol(nameof(sys), "_bernoulli")
    odes = ODESystem(clean_eqs, iv, vars, ps; name=odename)

    new_system = typeof(sys)(odes, field_values[2:end]...)

    new_system

end
