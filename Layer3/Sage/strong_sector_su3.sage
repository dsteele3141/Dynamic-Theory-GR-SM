#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE: Report (LA-8370-MS) Eqs 105-112 (5D Dirac, Williams' 8x8 rep, Eq
#   109 A-matrix corrected); master chart route 3 "charge quantisation";
#   book (2011) Eqs 7.66-7.67 (q = Z e, Z = N0 - N1 an integer) and the
#   non-singular potential phi_0 = k e^{-lambda/r}/r.  Williams' own proton
#   (two positrons + one electron, book Ch. 9) is deliberately NOT used.
# Context: weak_chirality_u4 (4x4 vs 8x8), weak_neutral_current (ledger Q6),
#   mass_term_su2 (Q12/Q13), hydrogen_nonsingular_spectrum (Q7 lambda bound).
# ---------------------------------------------------------------------------
#
# strong_sector_su3.sage                               [A] (assertion-hardened)
#
# PURPOSE.  Which parts of SU(3)_C can be DERIVED from DT plus the electroweak
# structure the data already forced on it (Q6, Q12, Q13), and which must be
# imported?
#
#   [A] room in DT's spinor: internal symmetry of the 5D Dirac field
#       (commutant of the Clifford algebra) -- 4x4: U(1); Williams' 8x8:
#       U(1)xU(1); an SU(N) needs N copies of ONE irrep
#   [B] consistency: DT's leptons + imported SU(2)_L x U(1)_Y are ANOMALOUS;
#       a doublet-pattern completion needs N_c copies with charges
#       Q_u = 1/2 + 1/(2 N_c), Q_d = -1/2 + 1/(2 N_c), and N_c ODD (Witten)
#   [C] DT's charge quantisation (q = Z e) then forbids a FREE quark: only
#       integer-charge composites can be isolated
#   [D] data fix N_c = 3: nucleon charges from three constituents, R ratio
#       below charm, Z hadronic width R_l (pi0 -> gamma gamma cannot: N_c drops out)
#   [E] dynamics: DT's own short-range potential phi_0 as the colour force --
#       running coupling and confinement against data
#
# Measured (imported, rounded): alpha_s(m_tau) = 0.330, alpha_s(M_Z) = 0.1180,
# alpha_s(~1 TeV) ~ 0.088 (PDG / CMS jets); R(2-3 GeV) = 2.2(1) (BES, KEDR);
# R_l = Gamma_had/Gamma_ll = 20.767(25) (LEP); GLS sum rule ~ 3 valence
# (CCFR); string tension sqrt(sigma) ~ 0.42 GeV (quarkonium, lattice); no free
# fractional charge in Millikan-type searches (Perl et al. 2009 review).
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("strong_sector_su3.sage -- what of SU(3)_C does DT supply?")
print("=" * 74)

K = CyclotomicField(4); i_ = K.gen()
def Mx(rows): return matrix(K, rows)
I2, Z2 = identity_matrix(K, 2), zero_matrix(K, 2)
s1, s2, s3 = Mx([[0, 1], [1, 0]]), Mx([[0, -i_], [i_, 0]]), Mx([[1, 0], [0, -1]])
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)

# ============================================================ [A]
print("\n[A] Internal symmetry room in DT's 5D Dirac field")
# 4x4 irrep (kappa = +1): Dirac gammas + gamma^4 = i gamma^5
g = [B(I2, Z2, Z2, -I2)] + [B(Z2, s, -s, Z2) for s in (s1, s2, s3)]
g5 = i_ * g[0] * g[1] * g[2] * g[3]
G4 = g + [i_ * g5]
# Williams' 8x8 (Eq 109 corrected; built exactly as in weak_chirality_u4)
Z4 = zero_matrix(K, 4)
beta = B(I2, Z2, Z2, -I2)
alph = [B(Z2, s, s, Z2) for s in (s1, s2, s3)]
Amat = B(Z2, i_ * I2, -i_ * I2, Z2)
bT = B(beta, Z4, Z4, beta)
aT = [B(Z4, a, a, Z4) for a in alph] + [B(Amat, Z4, Z4, Amat)]
G8 = [bT] + [-bT * a for a in aT]
eta = [1, -1, -1, -1, -1]
check(all(G8[j] * G8[k] + G8[k] * G8[j] == (2 * eta[j] if j == k else 0) * identity_matrix(K, 8)
          for j in range(5) for k in range(5)), "Williams' 8x8 satisfies the 5D Clifford algebra")
copies = lambda n: [block_diagonal_matrix([x] * n, subdivide=False) for x in G4]

def commutant(gens):
    n = gens[0].nrows(); In = identity_matrix(K, n)
    M = block_matrix([[x.transpose().tensor_product(In) - In.tensor_product(x)] for x in gens],
                     subdivide=False)
    basis = M.right_kernel().basis()
    return [matrix(K, n, n, list(v)).transpose() for v in basis]   # column-stacked vec
def abelian(C): return all(a * b == b * a for a in C for b in C)

rows = [("4x4 irrep (Q3, our choice)", G4, "U(1)"),
        ("Williams' 8x8 (Eq 109)", G8, "U(1) x U(1)"),
        ("2 copies of the 4x4", copies(2), "U(2) > SU(2)"),
        ("3 copies of the 4x4", copies(3), "U(3) > SU(3)")]
dims = {}
for name, gens, grp in rows:
    C = commutant(gens); dims[name] = (len(C), abelian(C))
    print("    %-28s commutant dim %d (%s)  ->  internal group %s"
          % (name, len(C), "abelian" if abelian(C) else "non-abelian", grp))
check(dims[rows[0][0]] == (1, True), "one 5D Dirac irrep has only a U(1) phase: no colour")
check(dims[rows[1][0]] == (2, True),
      "Williams' 8x8 = 4 + 4' (inequivalent): abelian U(1)xU(1), no SU(N) at all")
check(dims[rows[3][0]] == (9, False),
      "an SU(3) needs THREE copies of one irrep: colour is an ADDED multiplicity")
print("    DT's spinor (either rep) has no room for colour; like the SM, DT would")
print("    have to postulate the 3-fold multiplicity.")

# ============================================================ [B]
print("\n[B] Anomalies: DT's leptons + the imported SU(2)_L x U(1)_Y")
N, y = var('N y')
# Q = T3 + Y.  Higgs Y = 1/2 (Q12); Yukawa pattern (Q13) fixes
# Y_u = Y_Q + 1/2, Y_d = Y_Q - 1/2, Y_e = Y_L - 1/2.  Leptons: Q_e = -1 -> Y_L = -1/2.
YL = -QQ(1) / 2; Ye = YL - QQ(1) / 2
def anomalies(Nc, YQ):
    Yu, Yd = YQ + QQ(1) / 2, YQ - QQ(1) / 2
    su2 = Nc * YQ + YL                                         # SU(2)^2 U(1)
    y3 = Nc * (2 * YQ**3 - Yu**3 - Yd**3) + 2 * YL**3 - Ye**3  # U(1)^3 (L minus R)
    grav = Nc * (2 * YQ - Yu - Yd) + 2 * YL - Ye               # grav^2 U(1)
    return su2, y3, grav
a0 = anomalies(0, 0)
print("    leptons only (nu_L zero mode + e): SU(2)^2Y = %s, Y^3 = %s, grav = %s" % a0)
check(a0[0] != 0 and a0[1] != 0,
      "DT's leptons alone + imported SU(2)_L x U(1)_Y are ANOMALOUS: fermions must be added")
sol = solve(anomalies(N, y)[0] == 0, y)[0].rhs()
print("    add N_c copies of a doublet + two singlets: SU(2)^2Y = 0 forces Y_Q = %s" % sol)
check(bool(simplify(anomalies(N, sol)[1]) == 0) and bool(simplify(anomalies(N, sol)[2]) == 0),
      "the Y^3 and gravitational anomalies then cancel for EVERY N_c")
Qu, Qd = (QQ(1) / 2 + sol), (-QQ(1) / 2 + sol)
print("    charges: Q_u = %s, Q_d = %s" % (Qu.simplify_full(), Qd.simplify_full()))
check(all(bool(Qu.subs(N=n) != floor(Qu.subs(N=n))) for n in (3, 5, 7)),
      "for N_c > 1 the added fermions carry FRACTIONAL charge")
witten_ok = lambda n: (1 + n) % 2 == 0          # even number of SU(2) doublets
print("    Witten SU(2) anomaly: 1 lepton + N_c quark doublets must be even")
check([n for n in range(1, 8) if witten_ok(n)] == [1, 3, 5, 7], "=> N_c is ODD")

# ============================================================ [C]
print("\n[C] DT's charge quantisation (q = Z e, Eqs 7.66-7.67) on these charges")
def free_ok(nu, nd, n):          # net u, d counts (antiquark = -1); isolated => integer charge
    q = nu * Qu.subs(N=n) + nd * Qd.subs(N=n)
    return q == floor(q)
check(not free_ok(1, 0, 3) and not free_ok(0, 1, 3),
      "N_c = 3: a single quark (2/3 or -1/3) cannot be an isolated state")
check(free_ok(1, -1, 3) and free_ok(2, 1, 3) and free_ok(1, 2, 3) and not free_ok(1, 1, 3),
      "q-qbar mesons and 3-quark baryons are integer: only composites may be free")
print("    So DT's quantisation REQUIRES confinement of fractional charge (a")
print("    selection rule) -- it does not supply the mechanism; see [E].")

# ============================================================ [D]
print("\n[D] Data fix N_c")
for n in (1, 3, 5):
    qp = 2 * Qu.subs(N=n) + Qd.subs(N=n); qn = Qu.subs(N=n) + 2 * Qd.subs(N=n)
    print("    N_c = %d: Q_u = %s, Q_d = %s; three-constituent p, n charges = %s, %s"
          % (n, Qu.subs(N=n), Qd.subs(N=n), qp, qn))
check([n for n in range(1, 10, 2) if 2 * Qu.subs(N=n) + Qd.subs(N=n) == 1
       and Qu.subs(N=n) + 2 * Qd.subs(N=n) == 0] == [3],
      "nucleons = three valence constituents (GLS sum rule ~ 3) with charges +1, 0: only N_c = 3")
as_tau, as_Z = 0.330, 0.1180
Rm, Rs = 2.2, 0.1
Rf = lambda n: n * (Qu.subs(N=n)**2 + 2 * Qd.subs(N=n)**2) * (1 + as_tau / pi.n())
s2w = 0.2315
def Rl(n):
    gg = lambda T3, Q: (T3 - Q * s2w)**2 + (Q * s2w)**2
    up, dn = gg(0.5, Qu.subs(N=n)), gg(-0.5, Qd.subs(N=n))
    return n * (2 * up + 3 * dn) / gg(-0.5, -1) * (1 + as_Z / pi.n())
Rlm, Rls = 20.767, 0.025
for n in (1, 3, 5):
    print("    N_c = %d:  R(u,d,s) = %.2f (meas %.1f)   R_l(Z) = %.2f (meas %.3f)"
          % (n, Rf(n), Rm, Rl(n), Rlm))
check(abs(Rf(3) - Rm) < 2 * Rs and all(abs(Rf(n) - Rm) > 5 * Rs for n in (1, 5)),
      "R ratio below charm selects N_c = 3")
check(abs(Rl(3) - Rlm) / Rlm < 0.01 and all(abs(Rl(n) - Rlm) / Rlm > 0.2 for n in (1, 5)),
      "Z hadronic width selects N_c = 3 (tree x (1 + alpha_s/pi), < 1%)")
pi0 = [n * (Qu.subs(N=n)**2 - Qd.subs(N=n)**2) for n in (1, 3, 5, 7)]
check(all(bool(v == 1) for v in pi0),
      "pi0 -> gamma gamma: N_c (Q_u^2 - Q_d^2) = 1 for EVERY N_c once anomalies fix "
      "the charges -- it does not count colours (Baer-Wiese)")

# ============================================================ [E]
print("\n[E] DT's own short-range potential phi_0 as the colour force")
hbarc = 0.1973269804                               # GeV fm
# coupling at separation r: alpha_eff(r) = r^2 F(r)/hbar c ~ (1 - lam/r) e^{-lam/r}
f = lambda x: (1 - x) * exp(-x)                    # x = lam/r = lam Q/hbar c
m_tau, MZ, Q3 = 1.777, 91.1876, 1000.0
x_of = lambda lam, Q: lam * Q / hbarc
lam = find_root(lambda L: f(x_of(L, MZ)) / f(x_of(L, m_tau)) - as_Z / as_tau, 1e-6, 1.5e-3)
pred = as_Z * f(x_of(lam, Q3)) / f(x_of(lam, MZ))
print("    fit lambda to alpha(M_Z)/alpha(m_tau): lambda = %.2e fm" % lam)
print("    predicted alpha(1 TeV) = %+.4f   measured ~ +0.088" % pred)
check(pred < 0, "phi_0 running turns REPULSIVE by 1 TeV (core lam/r > 1): excluded")
print("    shape: phi_0 suppresses exponentially, e^{-lam Q}; QCD runs as 1/ln(Q^2/Lambda^2)")
r_big = 10.0                                         # fm
F_tail = (1 - lam / r_big) * exp(-lam / r_big) / r_big**2
print("    long range: phi_0 force ~ 1/r^2 -> separation energy FINITE (~ alpha hbar c / r)")
sigma = 0.42**2 / hbarc                              # GeV/fm, linear potential
check(F_tail * as_Z * hbarc < sigma / 100,
      "no string tension: phi_0 at 10 fm is < 1% of the measured confining force")
print("    => free fractional charges at finite energy, which [C] forbids:")
print("       DT's quantisation needs confinement that DT's potential cannot give.")

print("\n" + "=" * 74)
print("""VERDICT
  DERIVED (given the electroweak structure the data already forced, Q6/Q12/Q13):
  1. DT's leptons alone make the imported SU(2)_L x U(1)_Y anomalous ([B]);
     a doublet-pattern completion needs an N_c-fold multiplicity of fermions
     with charges Q_u = 1/2 + 1/(2N_c), Q_d = -1/2 + 1/(2N_c), and N_c odd.
     This is consistency (gauge anomalies), not DT physics -- the SM's own
     argument -- but it holds for DT once DT adopts the electroweak group.
  2. DT's charge quantisation (q = Z e) then forbids an isolated quark: only
     integer-charge composites can be free ([C]).  DT therefore REQUIRES
     confinement as a selection rule, and agrees with the null free-quark
     searches.
  SELECTED by data (not derived): N_c = 3 -- three-constituent nucleons with
  charges +1, 0, the R ratio, and the Z hadronic width all pick 3 ([D]).
  NOT SUPPLIED by DT:
  3. Room for colour: one 5D Dirac irrep has only a U(1) phase; Williams' 8x8
     is abelian U(1)xU(1).  SU(3) needs three copies of one irrep ([A]), so
     the multiplicity must be postulated -- as in the SM.
  4. The dynamics: DT's phi_0, used as the colour force, turns repulsive
     before 1 TeV (measured alpha_s ~ 0.09) and has no string tension ([E]).
     Gluons, the SU(3) gauge group and asymptotic freedom are imported.
  STANDING (assumption ledger: SHARED WITH SM, by import, with two derived
  constraints).  Just as the SM assumes quarks in the triplet of a gauged
  SU(3)_C, if DT assumes the same, its hadronic sector is consistent with
  the SM.  Ledger Q14 (colour multiplicity, N_c = 3: oddness and fractional
  charges forced by consistency, the value selected, not derived) and Q15
  (SU(3)_C gauge dynamics: imported; DT's phi_0 cannot serve).
  Bracket: (i) [B] assumes the completion follows the doublet + two singlets
  Yukawa pattern; other anomaly-free completions exist (e.g. vector-like
  pairs), so 'quarks' are the minimal choice, not the only one; (ii) [C]
  takes Williams' unit e (the electron's charge) as the quantum of Z e --
  with a unit e/3 free quarks would be allowed and [C] is lost; (iii) [E]
  applies Williams' NUCLEON potential to a role he never gave it (he does not
  use colour), so it shows DT supplies no colour force, NOT that his nuclear
  claim fails; Williams' e+e+e- proton is a separate, DT-specific test not
  run here; (iv) [D] uses imported QCD corrections (1 + alpha_s/pi), ~4-10%.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
