function construct_iter_all(N::Int, order::Int)

    #Construct an ordered iterator going over all moments
    # sequentially in terms of order

    # helper iterator through all possible combinations of
    # N terms each given in range [0, order]
    iter_full = Iterators.product(fill(0:order, N)...)
    iter_all = Iterators.filter(x -> sum(x) == 0, iter_full)
    for i in 1:order
        iter = Iterators.filter(x -> sum(x) == i, iter_full)
        iter_all = Iterators.flatten((iter_all, iter))
    end

    collect(iter_all)

end


# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))
# Function whichs run SymbolicUtils.simplify until there
# are no changes in the resulting expression (Operation)
#fully_simplify = Fixpoint(simplify);
# SymbolicUtils rule that expands a bracket
expansion_rule = @acrule ~x * +(~~ys) => sum(map(y-> simplify(~x * y), ~~ys));
# SymbolicUtils rule that expands all brackets and simplifies the expression until no changes occur
expand = Fixpoint(Prewalk(simplify(PassThrough(expansion_rule))));


PLUS_RULES = [
    @rule(~x::isnotflat(+) => flatten_term(+, ~x))
    @rule(~x::needs_sorting(+) => sort_args(+, ~x))
    @ordered_acrule(~a::is_literal_number + ~b::is_literal_number => ~a + ~b)

    @acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
    @acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...))
    @acrule(*(~~x, ~α) + *(~~x, ~β) => *(~α + ~β, (~~x)...))

    @acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
    @acrule(*(~α::is_literal_number, ~x) + ~x => *(~α + 1, ~x))
    @rule(+(~~x::hasrepeats) => +(merge_repeats(*, ~~x)...))

    @ordered_acrule((~z::_iszero + ~x) => ~x)
    @rule(+(~x) => ~x)
]

TIMES_RULES = [
    @rule(~x::isnotflat(*) => flatten_term(*, ~x))
    @rule(~x::needs_sorting(*) => sort_args(*, ~x))

    @ordered_acrule(~a::is_literal_number * ~b::is_literal_number => ~a * ~b)
    #@rule(*(~~x::hasrepeats) => *(merge_repeats(^, ~~x)...))

    @acrule((~y)^(~n) * ~y => (~y)^(~n+1))
    @ordered_acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

    @ordered_acrule((~z::_isone  * ~x) => ~x)
    @ordered_acrule((~z::_iszero *  ~x) => ~z)
    @rule(*(~x) => ~x)
]

ASSORTED_RULES = [
    @rule(identity(~x) => ~x)
    @rule(-(~x) => -1*~x)
    @rule(-(~x, ~y) => ~x + -1(~y))
    @rule(~x::_isone \ ~y => ~y)
    @rule(~x \ ~y => ~y / (~x))
    @rule(~x / ~y => ~x * pow(~y, -1))
    @rule(one(~x) => one(symtype(~x)))
    @rule(zero(~x) => zero(symtype(~x)))
    @rule(ifelse(~x::is_literal_number, ~y, ~z) => ~x ? ~y : ~z)
]

function modified_simplifier()
    rule_tree = [SymbolicUtils.If(istree, Chain(ASSORTED_RULES)),
                 SymbolicUtils.If(SymbolicUtils.is_operation(+),
                    Chain(PLUS_RULES)),
                 SymbolicUtils.If(SymbolicUtils.is_operation(*),
                    Chain(TIMES_RULES))] |> RestartedChain
    rule_tree
end

#simplify_mod(expr) = simplify(expr, rewriter=Fixpoint(Postwalk(modified_simplifier())))
simplify_mod = Fixpoint(Postwalk(PassThrough(modified_simplifier()))) # TODO: check that this modification still works for gamma closure
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> simplify_mod(~x * y), ~~ys));
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod))); #TODO: check if this is better than expand() for our needs
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x))
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod))
clean_expr = Fixpoint(Chain([simplify_mod, expand_mod, flatten_mod]))


# I the current implementation of ModelingToolkit.jl,
# the general symbolic type Num does not work well together with
# simplify() and other methods from SymbolicUtils...
# Although within the package's code there are certain methods
# alleviating the problem, here we use ModelingToolkit macros to initialize the
# variables as Num types (easier to define variables with subscripts) and then
# change them to SymbolicUtils.jl types such as Term{Real} .val of Num variable gives Term{Real}

function define_μ(N::Int, order::Int)

    iter = construct_iter_all(N, order)

    indices = []
    for i in iter
        push!(indices, trim_key(i))
    end

    @parameters t
    @variables μ[indices](t)

    μs = Dict()
    for (i, idx) in enumerate(iter)
        μs[idx] = μ[i].val
    end
    μ = μs

    μ[Tuple(zeros(Int, N))] = 1

    μ

end


function define_M(N::Int, order::Int)

    iter = construct_iter_all(N, order)

    indices = []
    for i in iter
        push!(indices, trim_key(i))
    end

    @parameters t
    @variables M[indices](t)

    Ms = Dict()
    for (i, idx) in enumerate(iter)
        if sum(idx)==0
            Ms[idx] = 1
        elseif sum(idx)==1
            Ms[idx] = 0
        else
            Ms[idx] = M[i].val
        end
    end
    M = Ms

    M

end
