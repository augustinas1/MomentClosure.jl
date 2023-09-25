using MomentClosure
using Test
using Latexify 
using Catalyst
using Distributions

@parameters p
@register_symbolic Distributions.Geometric(p)
m = rand(Distributions.Geometric(p)) # b - mean burst size => p = 1/(1+b)

rn = @reaction_network begin
      k_on*(1-g), 0 --> g
      k_off*P^2, g --> 0
      k_p, g --> g + $m*P
      γ_p, P --> 0
end
binary_vars = [1]

raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaws=false)
clean_eqs = bernoulli_moment_eqs(raw_eqs, binary_vars)
closed_raw_eqs = moment_closure(raw_eqs, "conditional gaussian", binary_vars)

# latexify output is very sensitive to Julia version and latexify + Symbolics updates...
expr = replace(raw"\begin{align*}
\frac{d\mu_{1 0}}{dt} =& k_{on}", "\r\n"=>"\n")
@test latexify(clean_eqs)[1:46] == expr

expr1 = replace(raw"\begin{align*}
\mu_{1 2} =& \frac{\mu_{1 1}^{2}}{\mu_{1 0}} \\
\mu_{1 3} =& \frac{3 \mu_{1 2} \mu_{1 1}}{\mu_{1 0}} + \frac{-2 \mu_{1 1}^{3}}{\mu_{1 0}^{2}}
\end{align*}
", "\r\n"=>"\n")
expr2 = replace(raw"\begin{align*}
\mu_{1 2} =& \frac{\mu_{1 1}^{2}}{\mu_{1 0}} \\
\mu_{1 3} =& \frac{-2 \mu_{1 1}^{3}}{\mu_{1 0}^{2}} + \frac{3 \mu_{1 1} \mu_{1 2}}{\mu_{1 0}}
\end{align*}
", "\r\n"=>"\n")
exprl = latexify(closed_raw_eqs, :closure)
@test (exprl == expr1) || (exprl == expr2)

@test_throws MethodError latexify(raw_eqs, :closure)