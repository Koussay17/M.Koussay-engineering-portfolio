#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Retraction capillaire et viscoelasticite.
(A) Loi de Taylor-Culick : solution exacte et le facteur sqrt(2) manquant
    dans le bilan d'energie.
(B) Oldroyd-B en extension uniaxiale : etat stationnaire exact,
    viscosite extensionnelle, et divergence a Wi = 1/2.
(C) FENE-P : regularisation par extensibilite finie.
(D) Formulation log-conforme : preservation de la positivite.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==================================================================== (A)
print("="*74)
print("(A) LOI DE TAYLOR-CULICK")
print("="*74)
sig, rho, h = 0.072, 1000.0, 1e-5          # eau, film de 10 microns
Vtc = np.sqrt(2*sig/(rho*h))
Ven = np.sqrt(4*sig/(rho*h))
print(f"  sigma={sig} N/m, rho={rho} kg/m3, h={h*1e6:.0f} microns")
print(f"  bilan de QUANTITE DE MOUVEMENT : V = sqrt(2 sigma/(rho h)) = {Vtc:.4f} m/s")
print(f"  bilan d'ENERGIE (faux)         : V = sqrt(4 sigma/(rho h)) = {Ven:.4f} m/s")
print(f"  rapport = {Ven/Vtc:.6f}  (exactement sqrt(2) = {np.sqrt(2):.6f})")
print("  -> la moitie de l'energie de surface est dissipee dans le bourrelet,")
print("     meme en fluide parfait : le raccord film/bourrelet est un choc.")

# transitoire : d(m V)/dt = 2 sigma, m = rho h X, V = dX/dt
def tc(t, y):
    X, mV = y
    V = mV/(rho*h*X) if X > 1e-12 else 0.0
    return [V, 2*sig]
s = solve_ivp(tc, [0, 2e-3], [1e-9, 0.0], rtol=1e-11, atol=1e-16,
              dense_output=True, method="DOP853")
tt = np.linspace(2e-5, 2e-3, 400)
X, mV = s.sol(tt); V = mV/(rho*h*X)
print(f"\n  integration transitoire : V(t) constante ?")
print(f"    V a t=0.1 ms : {np.interp(1e-4,tt,V):.6f} m/s")
print(f"    V a t=2.0 ms : {V[-1]:.6f} m/s")
print(f"    ecart a Taylor-Culick : {100*abs(V[-1]-Vtc)/Vtc:.4f} %")
print(f"  solution exacte : X = t sqrt(2 sigma/(rho h)) -> V constante des t=0.")

# ==================================================================== (B)
print()
print("="*74)
print("(B) OLDROYD-B EN EXTENSION UNIAXIALE")
print("="*74)
print("  L = diag(edot, -edot/2, -edot/2) ; A diagonal")
print("  dA_xx/dt = 2 edot A_xx - (A_xx-1)/lam")
print("  dA_yy/dt =  -edot A_yy - (A_yy-1)/lam")
print("\n  etat stationnaire exact : A_xx = 1/(1-2Wi), A_yy = 1/(1+Wi)")

lam = 1.0
def oldB(t, A, edot, lam):
    Axx, Ayy = A
    return [2*edot*Axx - (Axx-1)/lam, -edot*Ayy - (Ayy-1)/lam]

print(f"\n{'Wi':>7}{'A_xx num':>14}{'A_xx exact':>14}{'A_yy num':>12}{'A_yy exact':>12}")
for Wi in [0.1, 0.2, 0.3, 0.4, 0.45, 0.49]:
    edot = Wi/lam
    s = solve_ivp(oldB, [0, 400*lam], [1.0, 1.0], args=(edot, lam),
                  rtol=1e-12, atol=1e-14, method="Radau")
    Axx, Ayy = s.y[0,-1], s.y[1,-1]
    print(f"{Wi:7.2f}{Axx:14.6f}{1/(1-2*Wi):14.6f}{Ayy:12.6f}{1/(1+Wi):12.6f}")

print("\n  viscosite extensionnelle : eta_E = 3 eta_s + 3 eta_p/[(1-2Wi)(1+Wi)]")
etas, etap = 0.1, 0.9
print(f"{'Wi':>7}{'eta_E num':>14}{'eta_E exact':>14}")
for Wi in [0.1, 0.25, 0.4, 0.48]:
    edot = Wi/lam
    s = solve_ivp(oldB, [0, 400*lam], [1.0, 1.0], args=(edot, lam),
                  rtol=1e-12, atol=1e-14, method="Radau")
    Axx, Ayy = s.y[0,-1], s.y[1,-1]
    eE = 3*etas + (etap/lam)*(Axx-Ayy)/edot
    ex = 3*etas + 3*etap/((1-2*Wi)*(1+Wi))
    print(f"{Wi:7.2f}{eE:14.6f}{ex:14.6f}")

print("\n  au-dela du seuil : croissance exponentielle de taux (2Wi-1)/lam")
print(f"{'Wi':>7}{'taux mesure':>15}{'taux exact':>13}")
for Wi in [0.6, 0.8, 1.0, 1.5]:
    edot = Wi/lam
    s = solve_ivp(oldB, [0, 30*lam], [1.0, 1.0], args=(edot, lam),
                  rtol=1e-12, atol=1e-14, dense_output=True, method="Radau")
    t2 = np.linspace(20*lam, 30*lam, 200)
    A = s.sol(t2)[0]
    r = np.polyfit(t2, np.log(A), 1)[0]
    print(f"{Wi:7.2f}{r:15.8f}{(2*Wi-1)/lam:13.8f}")

# ==================================================================== (C)
print()
print("="*74)
print("(C) FENE-P : regularisation par extensibilite finie")
print("="*74)
print("  dA/dt = L A + A L^T - (f(A) A - I)/lam,  f = 1/(1 - tr A / Lmax^2)")
def fene(t, A, edot, lam, Lm2):
    Axx, Ayy = A
    tr = Axx + 2*Ayy
    f = 1.0/max(1 - tr/Lm2, 1e-12)
    return [2*edot*Axx - (f*Axx-1)/lam, -edot*Ayy - (f*Ayy-1)/lam]
print(f"{'Wi':>7}{'Lmax^2=100':>14}{'Lmax^2=1000':>14}{'Oldroyd-B':>14}")
for Wi in [0.3, 0.45, 0.5, 0.7, 1.0, 2.0]:
    edot = Wi/lam
    row = []
    for Lm2 in [100.0, 1000.0]:
        s = solve_ivp(fene, [0, 600*lam], [1.0,1.0], args=(edot,lam,Lm2),
                      rtol=1e-11, atol=1e-13, method="Radau")
        row.append(s.y[0,-1])
    ob = 1/(1-2*Wi) if Wi < 0.5 else np.inf
    obs = f"{ob:14.4f}" if np.isfinite(ob) else f"{'divergent':>14}"
    print(f"{Wi:7.2f}{row[0]:14.4f}{row[1]:14.4f}"+obs)
print("  -> FENE-P sature a A_xx < Lmax^2 : plus de singularite.")

# ==================================================================== (D)
print()
print("="*74)
print("(D) FORMULATION LOG-CONFORME : preservation de la positivite")
print("="*74)
print("  Cas critique : RELAXATION depuis un etat etire (apres un evenement")
print("  extensionnel), sans ecoulement.  dA/dt = -(A-1)/lam.")
print("  Euler explicite direct : A^(n+1) = A^n (1 - dt/lam) + dt/lam")
print("  -> si dt/lam > 1 et A^n grand, A^(n+1) devient NEGATIF (non physique).")
A0 = 50.0
print(f"\n  A_0 = {A0}, sans ecoulement")
print(f"{'dt/lam':>9}{'A_1 direct':>14}{'A>0 ?':>8}{'A_1 log-conf':>15}{'A>0 ?':>8}"
      f"{'nb pas avant A<0':>18}")
for r in [0.5, 0.9, 1.05, 1.5, 2.0, 4.0]:
    dt = r*lam
    A = A0; nneg = None
    for i in range(1, 2001):
        A = A + dt*(-(A-1)/lam)
        if i == 1: A1 = A
        if A <= 0 and nneg is None: nneg = i
        if not np.isfinite(A): break
    P = np.log(A0); ok = True
    for i in range(1, 2001):
        P = P + dt*(-(1-np.exp(-P))/lam)
        if i == 1: P1 = np.exp(P)
        if not np.isfinite(P): ok = False; break
    print(f"{r:9.2f}{A1:14.4f}{str(A1>0):>8}{P1:15.4f}{str(ok):>8}"
          f"{(str(nneg) if nneg else 'jamais'):>18}")
print("\n  -> A = exp(Psi) est positif par construction, quel que soit le pas.")
print("     C'est tout l'objet de la formulation de Fattal-Kupferman : la perte")
print("     de positivite du tenseur de conformation est LE mode de rupture des")
print("     simulations a grand nombre de Weissenberg.")

# ================================================= (E) nombre elastocapillaire
print()
print("="*74)
print("(E) QUAND LA SINGULARITE FRAPPE-T-ELLE LA RETRACTION ?")
print("="*74)
print("  taux d'elongation dans le bourrelet : edot ~ V/h avec V = Taylor-Culick")
print("  => Wi = lam V/h ;  la singularite Wi=1/2 est atteinte pour lam > h/(2V)")
for hh in [1e-6, 1e-5, 1e-4]:
    V = np.sqrt(2*sig/(rho*hh)); lc = hh/(2*V)
    print(f"    h={hh*1e6:6.1f} um : V={V:7.3f} m/s, edot={V/hh:.3e} 1/s, "
          f"lam_crit={lc:.3e} s")
print("  -> pour un film de 10 um, tout polymere de temps de relaxation")
print("     superieur a ~1 microseconde place la retraction au-dela du seuil :")
print("     Oldroyd-B y est structurellement inapplicable.")

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":.3,
                     "figure.dpi":150,"savefig.bbox":"tight","axes.linewidth":.8})
C1,C2,C3="#1f3b73","#b3452c","#2f7d4f"
fig,axs=plt.subplots(1,2,figsize=(9.4,3.3))
Wg=np.linspace(0.01,0.49,300)
axs[0].semilogy(Wg,3*etas+3*etap/((1-2*Wg)*(1+Wg)),color=C1,lw=1.8,
                label="Oldroyd-B (exact)")
for Lm2,c,ls in [(100.0,C2,"--"),(1000.0,C3,"-.")]:
    Wf=np.linspace(0.01,2.0,60); ee=[]
    for W in Wf:
        s=solve_ivp(fene,[0,600*lam],[1.,1.],args=(W/lam,lam,Lm2),
                    rtol=1e-10,atol=1e-12,method="Radau")
        Axx,Ayy=s.y[0,-1],s.y[1,-1]
        ee.append(3*etas+(etap/lam)*(Axx-Ayy)/(W/lam))
    axs[0].semilogy(Wf,ee,ls,color=c,lw=1.4,label=rf"FENE-P, $L^2={Lm2:.0f}$")
axs[0].axvline(0.5,color="k",ls=":",lw=1.2)
axs[0].annotate(r"$\mathrm{Wi}=1/2$",xy=(0.5,1e2),xytext=(0.56,3e2),fontsize=8)
axs[0].set_xlim(0,2.0); axs[0].set_xlabel(r"$\mathrm{Wi}=\lambda\dot\varepsilon$")
axs[0].set_ylabel(r"$\eta_E$"); axs[0].legend(fontsize=7.5)
axs[0].set_title("viscosite extensionnelle",fontsize=9.5)
for Wi2,c in zip([0.3,0.45,0.6,1.0],[C1,C3,C2,"#7a5aa8"]):
    s=solve_ivp(oldB,[0,25*lam],[1.,1.],args=(Wi2/lam,lam),
                rtol=1e-12,atol=1e-14,dense_output=True,method="Radau")
    t3=np.linspace(0,25*lam,400)
    axs[1].semilogy(t3/lam,s.sol(t3)[0],color=c,lw=1.5,label=rf"$\mathrm{{Wi}}={Wi2}$")
axs[1].set_xlabel(r"$t/\lambda$"); axs[1].set_ylabel(r"$A_{xx}$")
axs[1].legend(fontsize=8); axs[1].set_title(r"transitoire : borne si $\mathrm{Wi}<1/2$",fontsize=9.5)
fig.savefig("/home/claude/v_oldroyd.pdf")
print("\nfigure ecrite.")
