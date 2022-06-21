# determine ⟨xⁿ⟩, the n-th moment of the stoichiometric matrix coefficient x
# which is allowed to be a random variable
# NOTE: only independent geometrically distributed variables can be resolved
# TODO: follow changes in Catalyst and extend the compatibility accordingly
# TODO: extend for other univariate discrete distributions

expected_coeff(x, n::Int) = x^n # covers Number/Sym/Pow
expected_coeff(x::Mul, n::Int) = prod(expected_coeff(a, n) for a in arguments(x))

function expected_coeff(x::Div, n::Int)
    num, denom = arguments(x)
    (expected_coeff(num, n) / expected_coeff(denom, n))
end

expected_coeff(x::Add, n::Int) = sum( expected_coeff(a) for a in arguments(expand(x^n)) )
expected_coeff(x::Mul) = prod(expected_coeff(a) for a in arguments(x))
expected_coeff(x::Pow) = expected_coeff(arguments(x)...)
expected_coeff(x::Term) = expected_coeff(x, 1)
expected_coeff(x) = x

function expected_coeff(x::Term, n::Int)
    if isequal(operation(x), rand)
        args = arguments(x)
        isone(length(args)) || error("Unexpected arguments in $x")  
        resolve_moment(operation(args[1]), arguments(args[1]), n)
    else
        x^n
    end
end

function resolve_moment(d, args, n)
    @assert isequal(d, Geometric) "Can only compute moments for the geometric distribution."
    geometric_raw_moment(args[1], n)
end

eulerian_number(n::Int, k::Int) = sum( (-1)^j * binomial(n+1, j) * (k-j+1)^n for j in 0:k+1 )

function geometric_raw_moment(p, n::Int)
    # computes the symbolic expression of ⟨Xⁿ⟩ where X is a stochastic 
    # variable following a geometric distribution P(x) = (1-p)ˣp
    # ⟨Xⁿ⟩ = p*Li_{-k}(1-p) (polylogarithm of negative integer simplifies into a nicer sum)
    iszero(n) && return 1
    suma = 0.0
    for i in 0:n
        suma += eulerian_number(n, i) * (1-p)^(n-i)
        suma = expand(suma)
    end
    simplify(expand(p^(-n) * suma))
end