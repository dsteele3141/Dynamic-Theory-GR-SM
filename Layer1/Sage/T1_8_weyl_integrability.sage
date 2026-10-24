#!/usr/bin/env sage
# T1.8 - Weyl objection (second clock effect): two-sector escape
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), eq 65; T2.10.
# Supplementary Materials: step [L7], Sec. S2 (code index: Sec. S6).
#
# The Dynamic Theory is a Weyl theory (T1.8), so Einstein's 1918 objection
# applies prima facie: in Weyl geometry length transports as
#     d ell / ell = phi_i dx^i ,  ell(P) = ell0 exp( ∫ phi_i dx^i ),
# and a REAL phi with curl F_ij = d_i phi_j - d_j phi_i != 0 gives a path-
# dependent, history-dependent rescaling (atomic sizes/spectra would depend on
# worldline) — empirically excluded.
#
# This script machine-checks that the theory ESCAPES the objection in two
# complementary sectors (turning the escape from "argued" to "verified"):
#
#   A. GRAVITATIONAL sector is INTEGRABLE (Weyl-Integrable Spacetime). The gauge
#      vector is EXACT, phi_i = d_i (ln f^{1/2})  (eq 65), so its curl vanishes
#      identically and the closed-loop transport ∮ phi = 0 — no length
#      non-integrability, no second clock effect.
#
#   B. ELECTROMAGNETIC sector has an IMAGINARY coupling (alpha^1 = i alpha, T2.10),
#      so the Weyl factor is a U(1) PHASE, not a length rescaling: |ell| is
#      path-independent (only a phase accumulates) — the London-1927/Weyl-1929
#      reinterpretation. Real (observable) length is unchanged.
#
#   C. The accumulated phase is SINGLE-VALUED: the London quantisation
#      ∮ phi_k dx^k = 2*pi*i*n  gives  exp(∮ phi) = 1  (already proved in Lean,
#      T2_10_London.london_holonomy; re-confirmed here).
# ============================================================================

# Output is captured via stdout redirection by the runner.

_lines = []
def out(s=''):
    print(s); _lines.append(str(s))

out('=' * 70)
out('Weyl objection (second clock effect) — two-sector escape, verified')
out('=' * 70)
out('')

# ----------------------------------------------------------------------------
# A. Gravitational sector: exact gauge vector -> curl = 0 -> integrable (WIST)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(A) Gravitational sector INTEGRABLE: phi_i = d_i ln f^{1/2} exact => F=0')
out('-' * 70)
var('t x y z')
coords = [t, x, y, z]
f = function('f')(t, x, y, z)                 # arbitrary conformal factor (>0)
chi = (1/2) * log(f)                          # chi = ln f^{1/2}
phi = [diff(chi, ci) for ci in coords]        # phi_i = d_i chi   (eq 65, EXACT)

# Weyl curl F_ij = d_i phi_j - d_j phi_i must vanish identically (curl of grad)
allzero = True
for i in range(4):
    for j in range(4):
        Fij = (diff(phi[j], coords[i]) - diff(phi[i], coords[j])).simplify_full()
        if Fij != 0:
            allzero = False
            out('  F_%d%d = %s  [NONZERO]' % (i, j, Fij))
out('  F_ij = d_i phi_j - d_j phi_i = 0 for ALL i,j (any f): %s' % allzero)
assert allzero

# Closed-loop transport ∮ phi.dl = 0 (concrete: unit circle in (x,y), f>0)
out('  closed-loop check (concrete f = 1 + x^2 + y^2, unit circle):')
var('s')
fc = 1 + x^2 + y^2
chic = (1/2) * log(fc)
phix = diff(chic, x); phiy = diff(chic, y)
# parameterise x=cos s, y=sin s ; pull back the 1-form phi_x dx + phi_y dy
sub_loop = {x: cos(s), y: sin(s)}
integrand = (phix.subs(sub_loop) * (-sin(s)) + phiy.subs(sub_loop) * cos(s))
loop_int = integrate(integrand.simplify_full(), s, 0, 2*pi).simplify_full()
out('    ∮ phi_i dx^i = %s   (=0 ⇒ length integrable; no second clock effect)'
    % loop_int)
assert loop_int == 0
out('  => gravitational/conformal sector is Weyl-Integrable (WIST). [VERIFIED]')
out('')

# ----------------------------------------------------------------------------
# B. EM sector: imaginary coupling -> Weyl factor is a PHASE (|ell| invariant)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(B) EM sector: imaginary alpha^1 = i*alpha  =>  Weyl factor is a U(1) phase')
out('-' * 70)
var('a1 a2 a3 a4', domain='real')             # real EM gauge field components
phi_em = [I*a1, I*a2, I*a3, I*a4]             # phi imaginary (alpha^1 = i alpha)

# length: ell = ell0 exp(∫ phi.dx). Magnitude change d ln|ell| = Re(phi_i) dx^i.
re_parts = [real_part(p) for p in phi_em]
out('  Re(phi_i) = %s   (=> d ln|ell| = 0, magnitude path-independent)' % re_parts)
assert all(r == 0 for r in re_parts)

# the accumulated factor exp(i*theta) has modulus 1 (pure phase), theta real:
var('theta', domain='real')
zfac = exp(I*theta)
modsq = (zfac * zfac.conjugate()).simplify_full()      # |exp(i theta)|^2
out('  |exp(∫phi)|^2 = |exp(i*theta)|^2 = %s   (pure phase, no length rescaling)'
    % modsq)
assert modsq == 1
out('  => observable length |ell| is UNCHANGED; only an EM/U(1) phase transports.')
out('     (London 1927 / Weyl 1929 reinterpretation.)  [VERIFIED]')
out('')

# ----------------------------------------------------------------------------
# C. Holonomy single-valued: London quantisation  exp(∮phi)=exp(2*pi*i*n)=1
# ----------------------------------------------------------------------------
out('-' * 70)
out('(C) Single-valued phase: London quantisation exp(∮phi)=exp(2*pi*i*n)=1')
out('-' * 70)
ok_hol = all((exp(2*pi*I*nn)).simplify_full() == 1 for nn in [1, 2, 3, -1, -5])
out('  exp(2*pi*i*n) = 1 for n in {1,2,3,-1,-5}: %s' % ok_hol)
assert ok_hol
out('  (proved generally for integer n in Lean: T2_10_London.london_holonomy.)')
out('  => the EM phase is single-valued around closed loops: wavefunction well-')
out('     defined, no history dependence — exactly Bohr-Sommerfeld quantisation.')
out('')

# ----------------------------------------------------------------------------
out('=' * 70)
out('SUMMARY — the Dynamic Theory escapes Einstein\'s second-clock objection:')
out('  (A) gravity sector: exact phi (eq 65) => curl F = 0, ∮phi = 0 (WIST,')
out('      integrable) — no length non-integrability.                       [OK]')
out('  (B) EM sector: imaginary coupling (alpha^1=i alpha) => Weyl factor is a')
out('      pure U(1) phase, |ell| path-independent (London/Weyl-1929).       [OK]')
out('  (C) London quantisation => exp(∮phi)=1, phase single-valued.          [OK]')
out('')
out('  => the genuinely non-integrable part is the EM field, carried by an')
out('     imaginary coupling as an unobservable phase; the metric sector is')
out('     integrable. No REAL, observable second clock effect.')
out('')
out('  RESIDUAL CHECK (done in T1_8_conformal_gauge.sage): confirm the conformal')
out('  gauge ĝ -> f·ĝ introduces no real rescaling of physical lengths, and that')
out('  the exact piece (gravity) and non-exact F (EM) split is exactly Williams\'')
out('  decomposition. Reality structure here is consistent with that reading.')
out('=' * 70)

# (result captured via stdout redirection by the runner)
