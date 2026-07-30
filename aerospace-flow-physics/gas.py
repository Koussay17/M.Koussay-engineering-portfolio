#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lois de conservation hyperboliques et ecoulements compressibles.
(A) Burgers : catastrophe du gradient, et NON-UNICITE des solutions faibles.
(B) Condition d'entropie : les chocs de detente sont interdits (quantifie).
(C) Solveur de Riemann EXACT pour Euler : verification de Rankine-Hugoniot
    a la precision machine et production d'entropie.
(D) Tuyere quasi-1D : theoreme du col sonique, blocage, position du choc.
"""
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = 1.4                     # gamma
bet = (g-1)/(g+1)

# ==================================================================== (A)
print("="*74)
print("(A) BURGERS : CATASTROPHE DU GRADIENT ET NON-UNICITE")
print("="*74)
print("  u_t + u u_x = 0 ; u(x,0) = -sin(pi x) sur [0,2]")
print("  Les caracteristiques sont x = x0 + u0(x0) t : elles se croisent a")
print("  t_b = -1/min(u0') = 1/pi.")
tb = 1/np.pi
print(f"    t_b theorique = 1/pi = {tb:.8f}")
# detection numerique : jacobien 1 + t u0'(x0) = 0
x0 = np.linspace(0, 2, 200001)
u0p = -np.pi*np.cos(np.pi*x0)
tb_num = -1.0/u0p.min()
print(f"    t_b numerique (min du jacobien) = {tb_num:.8f}")
print(f"    ecart = {abs(tb_num-tb)/tb*100:.2e} %")

print("\n  NON-UNICITE : probleme de Riemann u_L=0, u_R=1 (profil croissant)")
uL, uR = 0.0, 1.0
print("    solution 1 : eventail de detente  u = x/t  pour 0<x<t")
print("    solution 2 : 'choc de detente' de vitesse s = (u_L+u_R)/2 = "
      f"{(uL+uR)/2:.3f}")
# les deux verifient Rankine-Hugoniot : s [u] = [u^2/2]
s2 = (uL+uR)/2
rh = s2*(uR-uL) - (uR**2/2 - uL**2/2)
print(f"    Rankine-Hugoniot pour la solution 2 : s[u]-[u^2/2] = {rh:.3e}  (=0)")
print("    -> les DEUX sont des solutions faibles. Rankine-Hugoniot ne suffit pas.")
print("    Condition de Lax (u_L > u_R pour un choc) : "
      f"{uL:.1f} > {uR:.1f} ? {uL>uR}  -> solution 2 REJETEE.")

# ==================================================================== (B)
print()
print("="*74)
print("(B) CONDITION D'ENTROPIE : les chocs de detente sont interdits")
print("="*74)
print("  Courbe de Hugoniot : rho2/rho1 = (P+beta)/(beta P+1),  P = p2/p1,")
print(f"  beta = (g-1)/(g+1) = {bet:.6f}")
print("  Saut d'entropie : ds/cv = ln P + g ln[(beta P + 1)/(P + beta)]")
def ds_cv(P):
    return np.log(P) + g*np.log((bet*P+1)/(P+bet))
print(f"\n{'P = p2/p1':>12}{'rho2/rho1':>12}{'ds/cv':>14}   verdict")
for P in [0.2, 0.5, 0.8, 0.95, 1.0, 1.05, 1.5, 3.0, 10.0]:
    r21 = (P+bet)/(bet*P+1)
    d = ds_cv(P)
    v = ("INTERDIT (ds<0)" if d < -1e-14 else
         "admissible" if d > 1e-14 else "limite (ds=0)")
    print(f"{P:12.2f}{r21:12.5f}{d:14.3e}   {v}")
print("\n  -> ds/cv > 0 si et seulement si P > 1 : seuls les chocs de")
print("     COMPRESSION sont admissibles. C'est le second principe qui")
print("     selectionne la solution parmi les solutions faibles.")
# comportement faiblement non lineaire
print("\n  regime faiblement non lineaire : le developpement de ds/cv en")
print("  P = 1+eps voit ses termes en eps et eps^2 s'annuler EXACTEMENT ;")
print("  il reste  ds/cv = (g^2-1)/(12 g^2) eps^3 + O(eps^4).")
cub = (g**2-1)/(12*g**2)
print(f"  coefficient exact = (g^2-1)/(12 g^2) = {cub:.8f}")
print(f"{'eps':>9}{'ds/cv exact':>16}{'cub*eps^3':>16}{'rapport':>10}")
for eps in [0.1, 0.03, 0.01, 0.003, 0.001]:
    d = ds_cv(1+eps)
    print(f"{eps:9.3f}{d:16.6e}{cub*eps**3:16.6e}{d/(cub*eps**3):10.5f}")
print("  -> un choc faible est presque isentropique : ds ~ eps^3.")

# ==================================================================== (C)
print()
print("="*74)
print("(C) SOLVEUR DE RIEMANN EXACT POUR EULER (tube a choc de Sod)")
print("="*74)
rL, uLv, pL = 1.0, 0.0, 1.0
rR, uRv, pR = 0.125, 0.0, 0.1
cL, cR = np.sqrt(g*pL/rL), np.sqrt(g*pR/rR)
print(f"  gauche : rho={rL}, u={uLv}, p={pL}   droite : rho={rR}, u={uRv}, p={pR}")

def fK(ps, rK, pK, cK):
    if ps > pK:                                   # choc
        A = 2/((g+1)*rK); B = bet*pK
        return (ps-pK)*np.sqrt(A/(ps+B)), np.sqrt(A/(ps+B))*(1-(ps-pK)/(2*(ps+B)))
    else:                                         # detente
        return (2*cK/(g-1))*((ps/pK)**((g-1)/(2*g))-1), (1/(rK*cK))*(ps/pK)**(-(g+1)/(2*g))

def phi(ps):
    a,_ = fK(ps, rL, pL, cL); b,_ = fK(ps, rR, pR, cR)
    return a + b + (uRv - uLv)

pstar = brentq(phi, 1e-8, 10.0, xtol=1e-15, rtol=1e-15)
aL,_ = fK(pstar, rL, pL, cL); aR,_ = fK(pstar, rR, pR, cR)
ustar = 0.5*(uLv+uRv) + 0.5*(aR-aL)
print(f"  p* = {pstar:.12f}   u* = {ustar:.12f}")
print(f"  (valeurs de reference de la litterature : 0.30313, 0.92745)")
# etats etoiles
rsL = rL*(pstar/pL)**(1/g)                       # detente a gauche
rsR = rR*((pstar/pR+bet)/(bet*pstar/pR+1))       # choc a droite
S = uRv + cR*np.sqrt((g+1)/(2*g)*pstar/pR + (g-1)/(2*g))
print(f"  rho*_L = {rsL:.10f} (detente)   rho*_R = {rsR:.10f} (choc)")
print(f"  vitesse du choc S = {S:.10f}")

print("\n  VERIFICATION DE RANKINE-HUGONIOT dans le repere du choc :")
v1, v2 = uRv - S, ustar - S
m1, m2 = rR*v1, rsR*v2
i1, i2 = pR + rR*v1**2, pstar + rsR*v2**2
h1 = g*pR/((g-1)*rR) + v1**2/2
h2 = g*pstar/((g-1)*rsR) + v2**2/2
print(f"    masse      : rho1 v1 = {m1:.12f}   rho2 v2 = {m2:.12f}   "
      f"ecart {abs(m1-m2):.2e}")
print(f"    impulsion  : p1+rho1 v1^2 = {i1:.12f}   p2+rho2 v2^2 = {i2:.12f}   "
      f"ecart {abs(i1-i2):.2e}")
print(f"    enthalpie  : h1+v1^2/2 = {h1:.12f}   h2+v2^2/2 = {h2:.12f}   "
      f"ecart {abs(h1-h2):.2e}")
dsR = np.log(pstar/pR) - g*np.log(rsR/rR)
dsL = np.log(pstar/pL) - g*np.log(rsL/rL)
print(f"\n    production d'entropie au choc  : ds/cv = {dsR:+.8f}  (> 0)")
print(f"    a travers la detente (gauche)  : ds/cv = {dsL:+.2e}  (= 0, isentropique)")
print(f"    -> le choc produit de l'entropie, la detente n'en produit pas.")

# ==================================================================== (D)
print()
print("="*74)
print("(D) TUYERE QUASI-1D")
print("="*74)
print("  dA/A = (M^2-1)/[M(1+(g-1)M^2/2)] dM")
print("  => M=1 exige dA=0 : le point sonique ne peut se trouver QU'AU COL.")
def AAs(M):
    return (1/M)*((2/(g+1))*(1+(g-1)/2*M**2))**((g+1)/(2*(g-1)))
print(f"\n  A/A* pour quelques M :")
for M in [0.2, 0.5, 0.9, 1.0, 1.1, 2.0, 3.0]:
    print(f"    M={M:4.1f} : A/A* = {AAs(M):.6f}")
print(f"  minimum de A/A* : {AAs(1.0):.10f} en M=1 (verification : "
      f"A/A*(0.999)={AAs(0.999):.8f}, A/A*(1.001)={AAs(1.001):.8f})")

print("\n  debit bloque : mdot sqrt(R T0)/(A* p0) = "
      f"{np.sqrt(g)*(2/(g+1))**((g+1)/(2*(g-1))):.8f}")

# choc dans une tuyere sur-detendue
Ae_At = 4.0
def M_super(AR):
    return brentq(lambda M: AAs(M)-AR, 1.0000001, 50.0)
def M_sub(AR):
    return brentq(lambda M: AAs(M)-AR, 1e-8, 0.9999999)
def p_p0(M): return (1+(g-1)/2*M**2)**(-g/(g-1))
Me_id = M_super(Ae_At)
print(f"\n  tuyere A_e/A_t = {Ae_At} : M_sortie ideal = {Me_id:.6f}, "
      f"p_e/p_0 = {p_p0(Me_id):.6f}")

def pe_for_shock(As_At):
    """Choc normal en A_s ; renvoie p_e/p_0."""
    M1 = M_super(As_At)
    M2 = np.sqrt((1+(g-1)/2*M1**2)/(g*M1**2-(g-1)/2))
    p2p1 = (2*g*M1**2-(g-1))/(g+1)
    p02p01 = p2p1*(1+(g-1)/2*M2**2)**(g/(g-1))/(1+(g-1)/2*M1**2)**(g/(g-1))
    # apres le choc : subsonique, A*_2 > A*_1
    As2_At = AAs(M2)                       # A_s / A*_2
    Ae_As2 = Ae_At/As_At*As2_At
    Me = M_sub(Ae_As2)
    return p_p0(Me)*p02p01, M1, M2, Me
print(f"\n{'A_s/A_t':>9}{'M1':>9}{'M2':>9}{'M_sortie':>11}{'p_e/p_0':>11}")
for As in [1.2, 1.6, 2.0, 2.8, 3.5, 3.99]:
    pe, M1, M2, Me = pe_for_shock(As)
    print(f"{As:9.2f}{M1:9.4f}{M2:9.4f}{Me:11.5f}{pe:11.6f}")
pb = 0.5
Asol = brentq(lambda A: pe_for_shock(A)[0]-pb, 1.0001, 3.999)
pe,M1,M2,Me = pe_for_shock(Asol)
print(f"\n  pour p_b/p_0 = {pb} : choc en A_s/A_t = {Asol:.6f} (M1={M1:.4f}), "
      f"verif p_e/p_0 = {pe:.8f}")

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":.3,
                     "figure.dpi":150,"savefig.bbox":"tight","axes.linewidth":.8})
C1,C2,C3="#1f3b73","#b3452c","#2f7d4f"
fig,axs=plt.subplots(1,3,figsize=(12,3.1))
Pg=np.logspace(np.log10(0.15),np.log10(12),400)
axs[0].semilogx(Pg,ds_cv(Pg),color=C1,lw=1.7)
axs[0].axhline(0,color="k",lw=.8); axs[0].axvline(1,color=C3,ls=":",lw=1.2)
axs[0].fill_between(Pg,ds_cv(Pg),0,where=ds_cv(Pg)<0,color=C2,alpha=.20)
axs[0].annotate("interdit",xy=(0.35,-0.02),fontsize=8,color=C2)
axs[0].annotate("admissible",xy=(3.0,0.05),fontsize=8,color=C1)
axs[0].set_xlabel(r"$P=p_2/p_1$"); axs[0].set_ylabel(r"$\Delta s/c_v$")
axs[0].set_title("condition d'entropie",fontsize=9.5)
# profil de Sod a t=0.2
t=0.2; x=np.linspace(-0.5,0.5,2000)
rho=np.empty_like(x); pp=np.empty_like(x); uu=np.empty_like(x)
csL=cL*(pstar/pL)**((g-1)/(2*g))
for i,xi in enumerate(x):
    xt=xi/t
    if xt < -cL: rho[i],uu[i],pp[i]=rL,uLv,pL
    elif xt < ustar-csL:
        uu[i]=2/(g+1)*(cL+xt); cc=cL-(g-1)/2*uu[i]
        rho[i]=rL*(cc/cL)**(2/(g-1)); pp[i]=pL*(cc/cL)**(2*g/(g-1))
    elif xt < ustar: rho[i],uu[i],pp[i]=rsL,ustar,pstar
    elif xt < S: rho[i],uu[i],pp[i]=rsR,ustar,pstar
    else: rho[i],uu[i],pp[i]=rR,uRv,pR
axs[1].plot(x,rho,color=C1,lw=1.5,label=r"$\rho$")
axs[1].plot(x,pp,color=C2,lw=1.5,label=r"$p$")
axs[1].plot(x,uu,color=C3,lw=1.5,label=r"$u$")
axs[1].set_xlabel(r"$x$"); axs[1].legend(fontsize=8)
axs[1].set_title(rf"Sod exact, $t={t}$",fontsize=9.5)
Ms=np.linspace(0.05,3.5,600)
axs[2].plot(Ms,[AAs(m) for m in Ms],color=C1,lw=1.7)
axs[2].plot([1],[1],"o",color=C2,ms=7,mfc="white",mew=1.6)
axs[2].annotate(r"$M=1$ : minimum de $A$",xy=(1,1),xytext=(1.3,2.2),fontsize=8,
                arrowprops=dict(arrowstyle="->",lw=.8))
axs[2].set_ylim(0.8,6); axs[2].set_xlabel(r"$M$"); axs[2].set_ylabel(r"$A/A^*$")
axs[2].set_title("relation section-Mach",fontsize=9.5)
fig.savefig("/home/claude/c_gas.pdf")
print("\nfigure ecrite.")
