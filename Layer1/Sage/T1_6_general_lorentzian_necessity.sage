# T1.6 - General Lorentzian necessity (the cross-term reduction)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), stability form 5.47.
# Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
from sage.all import *

# The Sage companion to the T1.6 signature-necessity Lean trio. Those modules
# prove that the BLOCK form Q = (q^0)^2 - h(x), with h positive-definite, has
# signature (1,n). This script closes the step they assume: that the GENERAL
# Williams stability form, which carries cross-terms 2 h_{0a} A dq^0 dq^a, is
# congruent to that block form WITHOUT loss of the signature conclusion -- i.e.
# positive-definiteness of the full stability Hessian survives the reduction.
# A symbolic proof of Schur positive-definiteness plus a multi-seed /
# multi-dimension robustness sweep; the assembled signature (1,n) is asserted.
#
# STRUCTURE
#   (1) Symbolic completion of the square (3+1D, full cross-terms): the general
#       Hessian quadratic form  =  f dY0^2 + S_{ab} dq^a dq^b  exactly, where
#       S = h - (1/f) b b^T is the Schur complement and dY0 the completed q^0.
#   (2) Symbolic Schur determinant identity  det(H) = f det(S).
#   (3) Symbolic PROOF that S is positive-definite when H is: setting dY0 = 0
#       (i.e. dq^0 = -(b.dq)/f) reduces the full form to exactly q^T S q, so a
#       value of the positive-definite H is exhibited at a nonzero vector for
#       every nonzero spatial q -- hence q^T S q > 0 for all q != 0.
#   (4) Numeric robustness sweep: many random positive-definite Hessians, dims
#       4 and 5; for each, assert f>0, S positive-definite, and that the
#       assembled spacetime form diag(+f) (+) (-S) has signature exactly (1,n).
#
# THE METRIC IDENTIFICATION (what the signature claim rests on, Rpt. 16-19).
#   The stability Hessian H is positive-definite (T1.5).  The absolute-velocity
#   arc-length constraint (T1.2) singles out the completed q^0-direction as the
#   TIMELIKE coordinate and fixes its sign +f>0; the remaining spatial Schur
#   block enters the line element with the OPPOSITE sign (-S).  Given that
#   identification, the work below shows the signature is forced to (1,n) for
#   the GENERAL form, not merely the aligned/diagonal one.  (This script
#   verifies the algebraic reduction + definiteness transfer; the physical sign
#   identification is the I-step of the Layer-1 reading.)
# ============================================================================


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 60)
print("T1.6  GENERAL LORENTZIAN NECESSITY (cross-term reduction)")
print("=" * 60)

# ----------------------------------------------------------------------------
# (1) Symbolic completion of the square (3+1D, full cross-terms)
# ----------------------------------------------------------------------------
print("\n1. Symbolic completion of the square (general 3+1D Hessian)")
var('f h01 h02 h03 h11 h22 h33 h12 h13 h23 dq0 dq1 dq2 dq3')

# General stability quadratic form: positive-definite Hessian H acting on
# (dq0, dq1, dq2, dq3), with a distinguished time component f = h00 and
# cross-terms 2 dq0 (h0a dqa).
cross_terms  = 2 * dq0 * (h01 * dq1 + h02 * dq2 + h03 * dq3)
spatial_pure = (h11 * dq1^2 + h22 * dq2^2 + h33 * dq3^2
                + 2 * h12 * dq1 * dq2 + 2 * h13 * dq1 * dq3 + 2 * h23 * dq2 * dq3)
Q_full = f * dq0^2 + cross_terms + spatial_pure

# Completed time coordinate  dY0 = dq0 + (b . dq)/f
dY0 = dq0 + (h01 * dq1 + h02 * dq2 + h03 * dq3) / f

# Schur complement  S = h_spatial - (1/f) b b^T
S11 = h11 - h01^2 / f;        S22 = h22 - h02^2 / f;        S33 = h33 - h03^2 / f
S12 = h12 - h01 * h02 / f;    S13 = h13 - h01 * h03 / f;    S23 = h23 - h02 * h03 / f
schur_spatial = (S11 * dq1^2 + S22 * dq2^2 + S33 * dq3^2
                 + 2 * S12 * dq1 * dq2 + 2 * S13 * dq1 * dq3 + 2 * S23 * dq2 * dq3)

completed_form = f * dY0^2 + schur_spatial
resid = (Q_full - completed_form).expand().simplify_full()
ok("general cross-term form == f*dY0^2 + Schur spatial form (congruence identity)",
   resid.is_zero())

# ----------------------------------------------------------------------------
# (2) Symbolic Schur determinant identity  det(H) = f det(S)
# ----------------------------------------------------------------------------
print("\n2. Schur determinant identity")
H = matrix([[f,   h01, h02, h03],
            [h01, h11, h12, h13],
            [h02, h12, h22, h23],
            [h03, h13, h23, h33]])
S = matrix([[S11, S12, S13],
            [S12, S22, S23],
            [S13, S23, S33]])
det_resid = (H.det() - f * S.det()).expand().simplify_full()
ok("det(H) = f * det(S) holds symbolically", det_resid.is_zero())

# ----------------------------------------------------------------------------
# (3) Symbolic PROOF: H positive-definite  =>  S positive-definite
#
#     The identity (1) gives  q^T H q = f*dY0^2 + q_s^T S q_s  for any vector
#     q = (dq0, q_s).  For an ARBITRARY nonzero spatial q_s, choose
#     dq0 = -(b . q_s)/f so that dY0 = 0; then q != 0 and
#         q_s^T S q_s = q^T H q > 0   (H positive-definite).
#     Hence S is positive-definite.  We verify the reduction step symbolically:
#     substituting dq0 = -(b.q_s)/f into the FULL form yields exactly the Schur
#     spatial form.
# ----------------------------------------------------------------------------
print("\n3. Symbolic positive-definiteness transfer (set dY0 = 0)")
dq0_star = -(h01 * dq1 + h02 * dq2 + h03 * dq3) / f
reduced = Q_full.subs({dq0: dq0_star}).expand().simplify_full()
ok("Q_full | (dY0=0) == Schur spatial form  =>  q_s^T S q_s is a value of H",
   (reduced - schur_spatial).expand().simplify_full().is_zero())
print("    therefore: H pos-def  =>  q_s^T S q_s > 0 for all q_s != 0  =>  S pos-def.")

# ----------------------------------------------------------------------------
# (4) Numeric robustness sweep (many random PD Hessians; dims 4 and 5)
#     Asserts, for each: f > 0, Schur S positive-definite, and the assembled
#     spacetime form  diag(+f) (+) (-S)  has signature exactly (1, n).
# ----------------------------------------------------------------------------
print("\n4. Numeric robustness sweep (random positive-definite Hessians)")
EPS = 1e-9
N_SEEDS = 200

def signature(mat):
    """(#positive, #negative, #zero) eigenvalues of a symmetric RDF matrix."""
    npos = nneg = nzero = 0
    for ev in mat.eigenvalues():
        x = float(ev)
        if x > EPS:      npos += 1
        elif x < -EPS:   nneg += 1
        else:            nzero += 1
    return npos, nneg, nzero

for dim in (4, 5):
    n_space = dim - 1
    n_pass = 0
    for seed in range(N_SEEDS):
        set_random_seed(seed)
        A = matrix(RDF, dim, dim, lambda i, j: random() - 0.5)
        Hn = A.transpose() * A + 2 * identity_matrix(RDF, dim)   # symmetric PD
        # sanity: Hn must itself be PD
        if not all(float(e) > EPS for e in Hn.eigenvalues()):
            continue
        fv = Hn[0, 0]
        b  = matrix(RDF, n_space, 1, [Hn[0, j] for j in range(1, dim)])
        M  = Hn.submatrix(1, 1, n_space, n_space)
        Sn = M - (1 / fv) * (b * b.transpose())                  # Schur complement
        # assembled spacetime form: +f time block (+) -S spatial block
        G = block_diagonal_matrix(matrix(RDF, [[fv]]), -Sn, subdivide=False)
        sig = signature(G)
        assert fv > EPS, "f must be > 0 (seed %d, dim %d)" % (seed, dim)
        assert all(float(e) > EPS for e in Sn.eigenvalues()), \
            "Schur complement not positive-definite (seed %d, dim %d)" % (seed, dim)
        assert sig == (1, n_space, 0), \
            "spacetime signature %s != (1,%d,0) (seed %d, dim %d)" % (sig, n_space, seed, dim)
        n_pass += 1
    ok("dim %d: all %d random PD Hessians give f>0, S pos-def, signature (1,%d)"
       % (dim, n_pass, n_space), n_pass == N_SEEDS)

# ----------------------------------------------------------------------------
print("\nConclusion")
print("  The general Williams stability form with h_0a cross-terms is congruent")
print("  to the block form f*dY0^2 + S (parts 1-2), and positive-definiteness of")
print("  the full Hessian transfers to the Schur complement S (part 3, symbolic;")
print("  part 4, numeric sweep).  Under the absolute-velocity identification of")
print("  the completed time direction (+f) against the spatial block (-S), the")
print("  signature is forced to (1,n) -- Lorentzian -- for the GENERAL form, not")
print("  only the aligned/diagonal specialisation checked in T1_6_stability_metric.")
print("=" * 60)
