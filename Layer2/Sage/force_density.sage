#!/usr/bin/env sage
# Force-density five-vector K = Div(T) = (1/c) F.J (Eq 70, p.102)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III.
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S2, Appendix A [C5].

from sage.all import *

def ok(label, cond):
    print(f"  {label}: {'PASS' if cond else 'FAIL'}")
    assert cond

E1, E2, E3 = var('E1 E2 E3', domain='real')
B1, B2, B3 = var('B1 B2 B3', domain='real')
V0, V1, V2, V3 = var('V0 V1 V2 V3', domain='real')
rho, J1, J2, J3, J4 = var('rho J1 J2 J3 J4', domain='real')
c = var('c', domain='real')
assume(c > 0)

print("="*70)
print("13_force_density_adjudicated.sage")
print("="*70)

print("\n(A) Bianchi identity cancellation (Div T ~ J F)")
# Coordinates: x0=ict, x1=x, x2=y
x0, x1, x2 = var('x0 x1 x2', domain='real')
coords = [x0, x1, x2]
F01 = function('F01')(*coords); F02 = function('F02')(*coords); F12 = function('F12')(*coords)
F3 = matrix(SR, [[0, F01, F02], [-F01, 0, F12], [-F02, -F12, 0]])
F_sq_scalar = sum(F3[m,n]^2 for m in range(3) for n in range(3))
term1 = sum(diff(F3[0,l], coords[k]) * F3[l,k] for l in range(3) for k in range(3))
term2 = 1/4 * diff(F_sq_scalar, x0)
# Bianchi: d0(F12) + d2(F01) + d1(F20) = 0 => d0(F12) = d1(F02) - d2(F01)
residual = (term1 + term2).subs(diff(F12, x0) == diff(F02, x1) - diff(F01, x2)).simplify_full()
ok("Bianchi identity reduces Div(T) to divergence-only terms", residual == 0)

print("\n(B) Force Component Verification (Page 102)")
F = matrix(SR, [
    [ 0,     I*E1,  I*E2,  I*E3,  I*V0],
    [-I*E1,  0,     B3,   -B2,    V1  ],
    [-I*E2, -B3,    0,     B1,    V2  ],
    [-I*E3,  B2,   -B1,   0,      V3  ],
    [-I*V0, -V1,   -V2,  -V3,    0   ]
])
J_cov = vector(SR, [I*c*rho, J1, J2, J3, J4])
K = [(-sum(J_cov[k] * F[k,j] for k in range(5)) / c).simplify_full() for j in range(5)]

ok("K0 matches corrected Entropy force (i/c)(J.E + J4*V0)", 
   (K[0] - (I/c)*(J1*E1 + J2*E2 + J3*E3 + J4*V0)).simplify_full() == 0)
ok("K1 matches corrected Lorentz + V-field force", 
   (K[1] - (rho*E1 + (J2*B3 - J3*B2)/c + J4*V1/c)).simplify_full() == 0)
ok("K4 matches Mass force rho*V0 - (J.V)/c",
   (K[4] - (rho*V0 - (J1*V1 + J2*V2 + J3*V3)/c)).simplify_full() == 0)

print("\n(C) GENERAL 5D divergence: Div(T) = F.J for ALL components, under dF=0")
# Part (A) checked only the j=0 Bianchi cancellation on a 3D subsystem. This
# section verifies the FULL 5D claim: the divergence of the energy-momentum
# tensor T_jk reduces to the 5-force K_j = (1/c)F_jk J^k, for every j, under the
# Bianchi identity dF=0.
cc = list(var('y0 y1 y2 y3 y4'))   # 5D coordinates (x0..x4)
eE1=function('eE1')(*cc); eE2=function('eE2')(*cc); eE3=function('eE3')(*cc)
bB1=function('bB1')(*cc); bB2=function('bB2')(*cc); bB3=function('bB3')(*cc)
vV0=function('vV0')(*cc); vV1=function('vV1')(*cc); vV2=function('vV2')(*cc); vV3=function('vV3')(*cc)
Ff = matrix(SR, 5, 5, [
    [ 0,     I*eE1,  I*eE2,  I*eE3,  I*vV0],
    [-I*eE1, 0,      bB3,   -bB2,    vV1 ],
    [-I*eE2,-bB3,    0,      bB1,    vV2 ],
    [-I*eE3, bB2,   -bB1,    0,      vV3 ],
    [-I*vV0,-vV1,   -vV2,   -vV3,    0   ]])

# T_jk = (F^2)_jk + 1/4 delta_jk (F_st F_st)
Fsq = sum(Ff[s,t]^2 for s in range(5) for t in range(5))
Tf  = Ff*Ff + (1/4)*Fsq*identity_matrix(SR, 5)
# source J_k = d_m F_km ; the 5-force (F.J)_j
Jf  = [sum(diff(Ff[k,m], cc[m]) for m in range(5)) for k in range(5)]

# Bianchi substitutions, keyed on atomic components (solving dF=0 per triple).
# Self-verification first: assert this set equals the genuine dF=0 3-form
# components auto-computed from F (so the reduction is neither over- nor
# under-strong) -- a full 10-of-10 match across all C(5,3) triples.
bianchi_subs = {
    diff(bB3,cc[0]): I*(diff(eE2,cc[1]) - diff(eE1,cc[2])),
    diff(bB2,cc[0]): I*(diff(eE1,cc[3]) - diff(eE3,cc[1])),
    diff(vV1,cc[0]): I*(diff(vV0,cc[1]) - diff(eE1,cc[4])),
    diff(bB1,cc[0]): I*(diff(eE3,cc[2]) - diff(eE2,cc[3])),
    diff(vV2,cc[0]): I*(diff(vV0,cc[2]) - diff(eE2,cc[4])),
    diff(vV3,cc[0]): I*(diff(vV0,cc[3]) - diff(eE3,cc[4])),
    diff(bB1,cc[1]): -diff(bB2,cc[2]) - diff(bB3,cc[3]),
    diff(vV2,cc[1]): diff(vV1,cc[2]) - diff(bB3,cc[4]),
    diff(vV3,cc[1]): diff(vV1,cc[3]) + diff(bB2,cc[4]),
    diff(vV3,cc[2]): diff(vV2,cc[3]) - diff(bB1,cc[4]),
}
auto_bianchi = []
for a in range(5):
    for b in range(a+1, 5):
        for d in range(b+1, 5):
            auto_bianchi.append((diff(Ff[b,d],cc[a]) + diff(Ff[d,a],cc[b])
                                 + diff(Ff[a,b],cc[d])).simplify_full())
n_matched = 0
for lhs, rhs in bianchi_subs.items():
    g = (lhs - rhs)
    if any(((g - B).simplify_full() == 0 or (g + B).simplify_full() == 0) for B in auto_bianchi):
        n_matched += 1
ok("the 10 Bianchi substitutions are exactly the dF=0 3-form components (10/10)",
   n_matched == 10 and len(bianchi_subs) == binomial(5, 3))

div_ok = True
for j in range(5):
    divT  = sum(diff(Tf[j,k], cc[k]) for k in range(5))
    force = sum(Ff[j,k]*Jf[k] for k in range(5))
    resid = (divT - force).subs(bianchi_subs).simplify_full()
    div_ok = div_ok and (resid == 0)
ok("Div(T)_j = (1/c-scaled) F_jk J^k for ALL j=0..4 under dF=0 (full 5D)", div_ok)

print("\nSummary: force-density components AND the full 5D divergence reduction")
print("Div(T)=F.J (all 5 components, under verified dF=0) machine-verified.")
print("="*70)
