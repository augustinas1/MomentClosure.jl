# determine the moment of the stoichiometric matrix coefficient
# necessary to deal with (independent) stochastic variables
# NOTE: only geometric distribution is supported
# the parameter passed is supposed to be the MEAN value
#function expected_coeff(x::Union{Int, Sym{ModelingToolkit.Parameter{Real}}}, pwr)
function expected_coeff(x::Union{Int, Symbolic}, pwr)
    if typeof(x) == Int64
        x^pwr
    else
        geometric_raw_moment(x, pwr)
    end
end

function eulerian_number(n, k)
    suma = 0
    for j in 0:k+1
        suma += (-1)^j * binomial(n+1, j) * (k-j+1)^n
    end
    return suma
end

#function geometric_raw_moment_arbitrary_order(m::Sym{ModelingToolkit.Parameter{Real}}, n::Int)
function geometric_raw_moment_arbitrary_order(m::Symbolic, n::Int)
    # computes the symbolic expression of ⟨Xⁿ⟩ where X is a stochastic variable
    # following a geometric distribution with mean burst size m
    # the distribution has the form P(x) = (1-p)ˣp, and m is given by m = (1-p)/p
    #@syms p
    suma = 0.0
    for i in 0:n
        #suma += eulerian_number(n, i) * (1-p)^(n-i)
        suma += eulerian_number(n, i) * m^(n-i) * (1+m)^(i-n)
        suma = simplify(suma)
    end
    #return p^(-n) * suma
    simplify((1+m)^n * suma)
end

#function geometric_raw_moment(m::Sym{ModelingToolkit.Parameter{Real}}, n::Int)
function geometric_raw_moment(m::Symbolic, n::Int)
    # provide the symbolic expression of ⟨Xⁿ⟩ where X is a stochastic variable
    # following a geometric distribution with mean burst size m
    # the distribution has the form P(x) = (1-p)ˣp, and m is given by m = (1-p)/p
    # Although geometric_raw_moment_arbitrary_order() can obtain the moment up to
    # arbitrary order, it is unable to format it in a nice and readable way
    # (not straightforward to use SymbolicUtils.jl to perform a binomial expansion
    # and hence tidy up the expression). For this reason, we provide the clean
    # symbolic expressions up to order 5 just to make the moment equations easier to check
    if n == 0
        1
    elseif n==1
        m
    elseif n==2
        2*m^2+m
    elseif n==3
        m*(6*m^2 + 6*m + 1)
    elseif n==4
        m*(2*m + 1)*(12*m^2+12*m+1)
    elseif n==5
        m*(120*m^4 + 240*m^3 + 150*m^2 + 30*m + 1)
    elseif n>5
        geometric_raw_moment_arbitrary_order(m, n)
    else
        error("negative argument values")
    end
end
