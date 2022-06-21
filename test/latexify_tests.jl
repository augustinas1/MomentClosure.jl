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
end k_on k_off k_p γ_p
binary_vars = [1]

raw_eqs = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=false)
clean_eqs = bernoulli_moment_eqs(raw_eqs, binary_vars)
closed_raw_eqs = moment_closure(raw_eqs, "conditional gaussian", binary_vars)

expr = replace(raw"\begin{align*}
\frac{d\mu_{1 0}}{dt} =& k_{on} - k_{on} \mu_{1 0} - k_{off} \mu_{1 2} \\
\frac{d\mu_{0 1}}{dt} =& k_{p} p^{-1} \mu_{1 0} - \gamma_{p} \mu_{0 1} - k_{p} \mu_{1 0} \\
\frac{d\mu_{1 1}}{dt} =& k_{on} \mu_{0 1} + k_{p} p^{-1} \mu_{1 0} - k_{p} \mu_{1 0} - k_{off} \mu_{1 3} - k_{on} \mu_{1 1} - \gamma_{p} \mu_{1 1} \\
\frac{d\mu_{0 2}}{dt} =& \gamma_{p} \mu_{0 1} + k_{p} \mu_{1 0} + 2 k_{p} p^{-2} \mu_{1 0} + 2 k_{p} p^{-1} \mu_{1 1} - 2 k_{p} \mu_{1 1} - 2 \gamma_{p} \mu_{0 2} - 3 k_{p} p^{-1} \mu_{1 0}
\end{align*}
", "\r\n"=>"\n")
@test latexify(clean_eqs) == expr

expr = replace(raw"\begin{align*}
\mu_{1 2} =& \mu_{1 0}^{-1} \mu_{1 1}^{2} \\
\mu_{1 3} =& 3 \mu_{1 0}^{-1} \mu_{1 1} \mu_{1 2} - 2 \mu_{1 0}^{-2} \mu_{1 1}^{3}
\end{align*}
", "\r\n"=>"\n")
@test latexify(closed_raw_eqs, :closure) == expr

@test_throws MethodError latexify(raw_eqs, :closure)