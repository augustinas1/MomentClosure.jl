using MomentClosure, Test, Latexify

@parameters t, k_on, k_off, k_p, γ_p, b
@variables p(t), g(t)

vars = [g, p]
ps = [k_on, k_off, k_p, γ_p, b]
S = [1 -1 0 0;
   0 0 b -1]
as = [k_on*(1-g),    # 0 -> g
      k_off*g*(p^2), # g -> 0
      k_p*g,         # 0 -> mP, m ~ Geometric(mean=b)
      γ_p*p]         # p -> 0
binary_vars = [1]
rn = ReactionSystemMod(t, vars, ps, as, S)

raw_eqs = generate_raw_moment_eqs(rn, 2)
closed_raw_eqs = moment_closure(raw_eqs, "conditional gaussian", binary_vars)

expr = replace(raw"\begin{align*}
\frac{d\mu_{1 0}}{dt} =& k_{on} - k_{on} \mu_{1 0} - k_{off} \mu_{1 0}^{-1} \mu_{1 1}^{2} \\
\frac{d\mu_{0 1}}{dt} =& b k_{p} \mu_{1 0} - \gamma_{p} \mu_{0 1} \\
\frac{d\mu_{1 1}}{dt} =& k_{on} \mu_{0 1} + b k_{p} \mu_{1 0} - k_{on} \mu_{1 1} - \gamma_{p} \mu_{1 1} - k_{off} \mu_{1 0}^{-2} \mu_{1 1}^{3} \\
\frac{d\mu_{0 2}}{dt} =& \gamma_{p} \mu_{0 1} + b k_{p} \mu_{1 0} + 2 k_{p} b^{2} \mu_{1 0} + 2 b k_{p} \mu_{1 1} - 2 \gamma_{p} \mu_{0 2}
\end{align*}
", "\r\n"=>"\n")
exprl = latexify(closed_raw_eqs)
@test exprl == expr

expr = replace(raw"\begin{align*}
\mu_{1 2} =& \mu_{1 0}^{-1} \mu_{1 1}^{2} \\
\mu_{1 3} =& 3 \mu_{1 0}^{-1} \mu_{1 1} \mu_{1 2} - 2 \mu_{1 0}^{-2} \mu_{1 1}^{3}
\end{align*}
", "\r\n"=>"\n")
@test latexify(closed_raw_eqs, :closure) == expr

@test_throws MethodError latexify(raw_eqs, :closure)