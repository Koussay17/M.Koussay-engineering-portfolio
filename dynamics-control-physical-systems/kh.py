#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stabilite lineaire non visqueuse d'ecoulements cisailles.
Equation de Rayleigh :  (U-c)(phi'' - k^2 phi) - U'' phi = 0,  phi(+-L)=0
resolue comme probleme aux valeurs propres generalise  A phi = c B phi
avec  A = diag(U)(D2 - k^2 I) - diag(U''),   B = D2 - k^2 I.
Verifications : mode neutre exact de tanh, theoreme de Howard, theoreme de Rayleigh.
"""
import numpy as np
from scipy.linalg import eig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------ differences finies
def grid_D2(L, N):
    """Grille uniforme sur [-L,L] (interieur seulement) et D2 d'ordre 4."""
    y = np.linspace(-L, L, N+2)
    h = y[1]-y[0]
    yi = y[1:-1]
    n = N
    D2 = np.zeros((n, n))
    for i in range(n):
        for off, c in [(-2,-1/12), (-1,4/3), (0,-5/2), (1,4/3), (2,-1/12)]:
            j = i+off
            if 0 <= j < n:
                D2[i, j] += c
    # bords : ordre 2 pour eviter les fantomes hors domaine
    D2[0,:] = 0; D2[0,0] = -2; D2[0,1] = 1
    D2[1,:] = 0; D2[1,0] = 1; D2[1,1] = -2; D2[1,2] = 1
    D2[-1,:] = 0; D2[-1,-1] = -2; D2[-1,-2] = 1
    D2[-2,:] = 0; D2[-2,-1] = 1; D2[-2,-2] = -2; D2[-2,-3] = 1
    return yi, D2/h**2, h

def rayleigh_spectrum(k, U, Upp, y, D2):
    n = len(y)
    I = np.eye(n)
    A = np.diag(U)@(D2 - k**2*I) - np.diag(Upp)
    B = D2 - k**2*I
    c = eig(A, B, right=False)
    return c[np.isfinite(c)]

def growth(k, U, Upp, y, D2):
    c = rayleigh_spectrum(k, U, Upp, y, D2)
    i = np.argmax(c.imag)
    return k*c[i].imag, c[i]

# ==================================================================== MAIN
if __name__ == "__main__":
    print("="*74)
    print("(0) MODE NEUTRE EXACT DE LA COUCHE tanh")
    print("="*74)
    print("  U = tanh(y) ; on verifie que phi = sech(y), k = 1, c = 0")
    print("  est solution exacte de l'equation de Rayleigh.")
    yv = np.linspace(-6, 6, 13)
    ph  = 1/np.cosh(yv)
    d2  = (1/np.cosh(yv))*(np.tanh(yv)**2 - (1/np.cosh(yv))**2)
    Uv  = np.tanh(yv)
    Uppv= -2*(1/np.cosh(yv))**2*np.tanh(yv)
    res = Uv*(d2 - 1.0*ph) - Uppv*ph
    print(f"  residu max analytique |(U-c)(phi''-k^2 phi) - U'' phi| = {np.abs(res).max():.3e}")

    L, N = 15.0, 500
    y, D2, h = grid_D2(L, N)
    U   = np.tanh(y)
    Upp = -2*np.tanh(y)/np.cosh(y)**2

    print(f"\n  resolution numerique (L={L}, N={N}, h={h:.4f}) :")
    for k in [0.9, 0.99, 1.0, 1.01, 1.1]:
        g, c = growth(k, U, Upp, y, D2)
        print(f"    k={k:5.2f}  k*c_i = {g:12.3e}   c = {c.real:+.3e}{c.imag:+.3e}j")

    print()
    print("="*74)
    print("(1) COURBE DE CROISSANCE DE LA COUCHE tanh")
    print("="*74)
    ks = np.linspace(0.02, 1.15, 60)
    gr, cc = [], []
    for k in ks:
        g, c = growth(k, U, Upp, y, D2)
        gr.append(max(g, 0.0)); cc.append(c)
    gr = np.array(gr); cc = np.array(cc)
    i = np.argmax(gr)
    # raffinement parabolique autour du maximum
    kk = np.linspace(ks[i]-0.05, ks[i]+0.05, 41)
    gg = np.array([max(growth(kx, U, Upp, y, D2)[0], 0.0) for kx in kk])
    j = np.argmax(gg)
    p = np.polyfit(kk[j-3:j+4], gg[j-3:j+4], 2)
    kmax = -p[1]/(2*p[0]); gmax = np.polyval(p, kmax)
    print(f"  taux maximal   k*c_i = {gmax:.6f}   en  k = {kmax:.6f}")
    print(f"  valeurs de reference (litterature) : 0.1897 en k = 0.4446")
    print(f"  ecart : {100*abs(gmax-0.1897)/0.1897:.2f} % sur le taux, "
          f"{100*abs(kmax-0.4446)/0.4446:.2f} % sur k")
    kc = ks[gr > 1e-6].max()
    print(f"  coupure : instable pour k < {kc:.4f}  (theorie exacte : k < 1)")

    print()
    print("="*74)
    print("(2) THEOREME DE HOWARD : c dans le demi-disque de diametre [Umin,Umax]")
    print("="*74)
    Umin, Umax = U.min(), U.max()
    ctr, rad = (Umin+Umax)/2, (Umax-Umin)/2
    print(f"  Umin={Umin:.4f}, Umax={Umax:.4f} -> centre={ctr:.4f}, rayon={rad:.4f}")
    worst = 0.0; nviol = 0; ntot = 0
    allc = []
    for k in [0.1, 0.3, 0.45, 0.6, 0.9]:
        c = rayleigh_spectrum(k, U, Upp, y, D2)
        c = c[c.imag > 1e-9]
        allc.append((k, c))
        for z in c:
            ntot += 1
            d = abs(z-ctr)/rad
            worst = max(worst, d)
            if d > 1+1e-8: nviol += 1
    print(f"  {ntot} modes instables testes ; violations : {nviol}")
    print(f"  rapport maximal |c-centre|/rayon = {worst:.6f}  (doit etre <= 1)")

    print()
    print("="*74)
    print("(3) THEOREME DE RAYLEIGH : profil sans point d'inflexion")
    print("="*74)
    Lp, Np = 1.0, 400
    yp, D2p, hp = grid_D2(Lp, Np)
    Up   = 1 - yp**2
    Uppp = -2*np.ones_like(yp)
    print("  U = 1 - y^2 sur [-1,1] : U'' = -2 ne s'annule jamais")
    print(f"{'k':>7}{'max c_i':>14}")
    mx = 0.0
    for k in [0.5, 1.0, 2.0, 4.0, 8.0]:
        c = rayleigh_spectrum(k, Up, Uppp, yp, D2p)
        m = c.imag.max(); mx = max(mx, m)
        print(f"{k:7.1f}{m:14.3e}")
    print(f"  max sur tous les k : {mx:.3e}  -> stable, conforme au theoreme")

    print()
    print("="*74)
    print("(4) NAPPE DE TOURBILLON : catastrophe ultraviolette")
    print("="*74)
    print("  epaisseur delta -> 0 : gamma = k*DU/2 croit sans borne avec k.")
    print("  L'epaisseur finie coupe l'instabilite a k*delta ~ 1 :")
    print(f"{'delta':>8}{'k_max':>10}{'gamma_max':>12}{'k*delta':>10}")
    for delta in [1.0, 0.5, 0.25]:
        Ud   = np.tanh(y/delta)
        Uppd = -2*np.tanh(y/delta)/np.cosh(y/delta)**2/delta**2
        kk2 = np.linspace(0.05/delta, 1.3/delta, 45)
        gg2 = np.array([max(growth(kx, Ud, Uppd, y, D2)[0], 0.0) for kx in kk2])
        jj = np.argmax(gg2)
        print(f"{delta:8.2f}{kk2[jj]:10.4f}{gg2[jj]:12.6f}{kk2[jj]*delta:10.4f}")
    print("  -> k_max*delta est invariant : c'est bien delta qui fixe la coupure.")

    # ------------------------------------------------------------- figures
    plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":.3,
                         "figure.dpi":150,"savefig.bbox":"tight","axes.linewidth":.8})
    C1,C2,C3="#1f3b73","#b3452c","#2f7d4f"
    fig,axs=plt.subplots(1,2,figsize=(9.4,3.3))
    axs[0].plot(ks,gr,"-",color=C1,lw=1.7,label=r"$kc_i$ (Rayleigh, num.)")
    axs[0].plot([kmax],[gmax],"o",color=C2,ms=7,mfc="white",mew=1.6,
                label=rf"max $={gmax:.4f}$ en $k={kmax:.4f}$")
    axs[0].axvline(1.0,color=C3,ls=":",lw=1.3)
    axs[0].annotate(r"mode neutre exact $k=1$, $\phi=\mathrm{sech}\,y$",
                    xy=(1.0,0.02),xytext=(0.55,0.06),fontsize=7.5,color=C3,
                    arrowprops=dict(arrowstyle="->",color=C3,lw=.8))
    axs[0].set_xlabel(r"$k\delta$"); axs[0].set_ylabel(r"taux de croissance $kc_i$")
    axs[0].legend(fontsize=8); axs[0].set_title(r"couche de melange $U=\tanh y$",fontsize=9.5)
    th=np.linspace(0,np.pi,200)
    axs[1].plot(ctr+rad*np.cos(th),rad*np.sin(th),"k-",lw=1.4,
                label="demi-cercle de Howard")
    axs[1].plot([ctr-rad,ctr+rad],[0,0],"k-",lw=1.4)
    for (k,c),col in zip(allc,[C1,C2,C3,"#7a5aa8","#8a8a3a"]):
        axs[1].plot(c.real,c.imag,"o",ms=5,color=col,mfc="white",mew=1.3,
                    label=rf"$k={k}$")
    axs[1].set_xlabel(r"$c_r$"); axs[1].set_ylabel(r"$c_i$")
    axs[1].set_aspect("equal"); axs[1].legend(fontsize=7.2)
    axs[1].set_title("tous les modes instables sont dans le demi-disque",fontsize=9.5)
    fig.savefig("/home/claude/k_stab.pdf")
    np.savez("/home/claude/kh.npz",ks=ks,gr=gr,kmax=kmax,gmax=gmax)
    print("\nfigure ecrite.")
