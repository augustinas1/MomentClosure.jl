using MomentClosure, Symbolics
using Test

@variables x, y, z
@variables u, v, w
@parameters p, q, r

ps = [p]
eq = p*x
subs = Dict(x => u)
sub_eq = MomentClosure.poly_subs(eq, subs, ps)

@test isequal(sub_eq, p*u)

ps = [p, q, r]
eq = ((p*x + q*y)^2 + r*z)^2
subs = Dict(y^4 => u, x*y*z => v, z*x^2 => w)
sub_eq = MomentClosure.poly_subs(eq, subs, ps)

@test isequal(sub_eq, (p^4)*(x^4) + (q^4)*u + 
                      (r^2)*(z^2) + 4p*x*(q^3)*(y^3) + 
                      4q*y*(p^3)*(x^3) + 2r*(p^2)*w + 2r*z*(q^2)*(y^2) + 
                      6(p^2)*(q^2)*(x^2)*(y^2) + 4p*q*r*v)


@variables x, y, z
@variables u, v, w

f = x * y
n = 2
@test isequal(0, MomentClosure.nth_differential(f, x, n))
@test isequal(u*v + v * (x - u) + u * (y - v) + (x - u)*(y - v), MomentClosure.taylor_expand(f, [x, y], [u, v], n))
f = exp(x)
n = 10
@test isequal(exp(x), MomentClosure.nth_differential(f, x, 10))
@test isequal(sum(x^i/factorial(i) for i in 0:n), MomentClosure.taylor_expand(f, [x], [0], n))

f = 1/x
n = 3
@test isequal((-1)^n*factorial(n)/x^(n+1), MomentClosure.nth_differential(f, x, n))
@test isequal(sum((-1)^i*(x-1)^i for i in 0:n), MomentClosure.taylor_expand(f, [x], [1], n))
