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
\frac{d\mu{_{10}}}{dt} =& k_{on} - k_{on} \mu{_{10}} - k_{off} \mu{_{10}}^{-1} \mu{_{11}}^{2} \\
\frac{d\mu{_{01}}}{dt} =& b k_{p} \mu{_{10}} - \gamma_{p} \mu{_{01}} \\
\frac{d\mu{_{11}}}{dt} =& k_{on} \mu{_{01}} + b k_{p} \mu{_{10}} - k_{on} \mu{_{11}} - \gamma_{p} \mu{_{11}} - k_{off} \mu{_{10}}^{-2} \mu{_{11}}^{3} \\
\frac{d\mu{_{02}}}{dt} =& \gamma_{p} \mu{_{01}} + b k_{p} \mu{_{10}} + 2 k_{p} b^{2} \mu{_{10}} + 2 b k_{p} \mu{_{11}} - 2 \gamma_{p} \mu{_{02}}
\end{align*}
", "\r\n"=>"\n")
@test latexify(closed_raw_eqs) == expr

expr = replace(raw"\begin{align*}
\mu{_{12}} =& \mu{_{10}}^{-1} \mu{_{11}}^{2} \\
\mu{_{13}} =& 3 \mu{_{10}}^{-1} \mu{_{11}} \mu{_{12}} - 2 \mu{_{10}}^{-2} \mu{_{11}}^{3}
\end{align*}
", "\r\n"=>"\n")
@test latexify(closed_raw_eqs, :closure) == expr

@test_throws MethodError latexify(raw_eqs, :closure)