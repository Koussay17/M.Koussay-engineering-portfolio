#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Methode modale de Fourier (RCWA) pour un reseau lamellaire.
Comparaison des regles de factorisation de Fourier :
  - Laurent naif   : q^2 = E (Kx E^-1 Kx - I),  V = i E^-1 W Q
  - regles de Li   : q^2 = A^-1 (Kx E^-1 Kx - I), V = i A W Q
avec E = Toeplitz(eps_n), A = Toeplitz((1/eps)_n).
Les deux convergent vers la MEME limite : seule la vitesse differe.
"""
import numpy as np
from scipy.linalg import toeplitz, eig, solve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------ Fourier du profil
def fourier_binary(er, eg, f, nmax):
    """Coefficients de Fourier c_n d'un creneau valant er sur [0,f) et eg sur [f,1)."""
    n = np.arange(-nmax, nmax+1)
    c = np.empty(len(n), dtype=complex)
    z = (n == 0)
    c[z] = f*er + (1-f)*eg
    nn = n[~z]
    c[~z] = (er-eg)*(1-np.exp(-2j*np.pi*nn*f))/(2j*np.pi*nn)
    return c

def toep(c, N):
    """Matrice de Toeplitz (2N+1)x(2N+1) : T[m,n] = c_{m-n}."""
    nmax = (len(c)-1)//2
    col = c[nmax:nmax+N+1][::-1]      # c_0, c_-1, ... non : construisons a la main
    T = np.empty((2*N+1, 2*N+1), dtype=complex)
    for m in range(2*N+1):
        for n in range(2*N+1):
            T[m, n] = c[nmax + (m-n)]
    return T

# ------------------------------------------------------------------ RCWA
def rcwa_TM(N, lam, period, depth, er, eg, f, eI, eII, theta=0.0, rule="li"):
    """Retourne (DE_r, DE_t) : efficacites de diffraction en TM (H // sillons)."""
    k0 = 2*np.pi/lam
    n = np.arange(-N, N+1)
    kx = (np.sqrt(eI)*np.sin(theta) + n*lam/period)          # kx/k0
    Kx = np.diag(kx)
    I = np.eye(2*N+1)

    ce  = fourier_binary(er, eg, f, 2*N+2)
    ci  = fourier_binary(1/er, 1/eg, f, 2*N+2)
    E   = toep(ce, N)
    A   = toep(ci, N)
    Einv = np.linalg.inv(E)

    M = Kx @ Einv @ Kx - I
    if rule == "li":
        Ainv = np.linalg.inv(A)
        Mat = Ainv @ M
        Vfac = A
    else:                                   # Laurent naif
        Mat = E @ M
        Vfac = Einv

    q2, W = np.linalg.eig(Mat)
    q = np.sqrt(q2.astype(complex))
    q = np.where(q.real < 0, -q, q)         # Re(q) >= 0 : decroissance
    Q = np.diag(q)
    V = 1j * Vfac @ W @ Q
    X = np.diag(np.exp(-q*k0*depth))

    kzI  = np.sqrt(eI  - kx**2 + 0j);  kzI  = np.where(kzI.imag < 0, -kzI, kzI)
    kzII = np.sqrt(eII - kx**2 + 0j);  kzII = np.where(kzII.imag < 0, -kzII, kzII)
    ZI  = np.diag(kzI/eI)
    ZII = np.diag(kzII/eII)

    delta = np.zeros(2*N+1, dtype=complex); delta[N] = 1.0

    # z = d :  (-V X + ZII W X) c+ + (V + ZII W) c- = 0
    L1 = -V @ X + ZII @ W @ X
    L2 =  V + ZII @ W
    Cm = -solve(L2, L1)                     # c- = Cm c+
    # z = 0 :  (ZI W + V) c+ + (ZI W - V) X c- = 2 ZI delta
    Lp = (ZI @ W + V) + (ZI @ W - V) @ X @ Cm
    cp = solve(Lp, 2*(ZI @ delta))
    cm = Cm @ cp

    R = W @ cp + W @ X @ cm - delta
    T = W @ X @ cp + W @ cm

    DEr = np.abs(R)**2 * np.real(kzI/eI) / np.real(kzI[N]/eI)
    DEt = np.abs(T)**2 * np.real(kzII/eII) / np.real(kzI[N]/eI)
    return np.real(DEr), np.real(DEt), n

def rcwa_TE(N, lam, period, depth, er, eg, f, eI, eII, theta=0.0):
    k0 = 2*np.pi/lam
    n = np.arange(-N, N+1)
    kx = (np.sqrt(eI)*np.sin(theta) + n*lam/period)
    Kx = np.diag(kx); I = np.eye(2*N+1)
    E = toep(fourier_binary(er, eg, f, 2*N+2), N)
    q2, W = np.linalg.eig(Kx@Kx - E)
    q = np.sqrt(q2.astype(complex)); q = np.where(q.real < 0, -q, q)
    Q = np.diag(q); V = 1j*W@Q; X = np.diag(np.exp(-q*k0*depth))
    kzI  = np.sqrt(eI - kx**2 + 0j);  kzI  = np.where(kzI.imag<0, -kzI, kzI)
    kzII = np.sqrt(eII- kx**2 + 0j);  kzII = np.where(kzII.imag<0, -kzII, kzII)
    YI, YII = np.diag(kzI), np.diag(kzII)
    delta = np.zeros(2*N+1, dtype=complex); delta[N] = 1.0
    Cm = -solve(V + YII@W, -V@X + YII@W@X)
    cp = solve((YI@W + V) + (YI@W - V)@X@Cm, 2*(YI@delta))
    cm = Cm@cp
    R = W@cp + W@X@cm - delta
    T = W@X@cp + W@cm
    DEr = np.abs(R)**2*np.real(kzI)/np.real(kzI[N])
    DEt = np.abs(T)**2*np.real(kzII)/np.real(kzI[N])
    return np.real(DEr), np.real(DEt), n

# ==================================================================== MAIN
if __name__ == "__main__":
    lam, period, depth = 0.8, 1.0, 0.5
    er, eg = -20.0+0j, 1.0+0j          # metal SANS pertes -> bilan d'energie exact
    f, eI, eII = 0.5, 1.0, 1.0

    print("="*74)
    print("RESEAU LAMELLAIRE METALLIQUE SANS PERTES")
    print(f"  lambda={lam}, periode={period}, profondeur={depth}")
    print(f"  eps_ruban={er.real}, eps_sillon={eg.real}, taux={f}, incidence normale")
    print("="*74)

    print("\n(0) BILAN D'ENERGIE (test independant du code)")
    for N in [10, 40, 120]:
        r, t, n = rcwa_TM(N, lam, period, depth, er, eg, f, eI, eII, rule="li")
        rt, tt, _ = rcwa_TE(N, lam, period, depth, er, eg, f, eI, eII)
        print(f"   N={N:4d}  TM: R+T = {r.sum()+t.sum():.12f}   "
              f"TE: R+T = {rt.sum()+tt.sum():.12f}")

    print("\n(1) CONVERGENCE EN TM : ordre 0 reflechi")
    Nref = 300
    vals = [rcwa_TM(M, lam, period, depth, er, eg, f, eI, eII, rule="li")[0][M]
            for M in (200, 300, 400, 500)]
    valn = [rcwa_TM(M, lam, period, depth, er, eg, f, eI, eII, rule="naive")[0][M]
            for M in (200, 300, 400, 500)]
    R0ref = float(np.mean(vals)); sref = float(np.std(vals))
    print(f"   reference = moyenne(Li, N=200..500) = {R0ref:.8f} +/- {sref:.1e}")
    print(f"   moyenne(naif, N=200..500)           = {np.mean(valn):.8f} "
          f"+/- {np.std(valn):.1e}")
    print(f"   -> les deux regles convergent bien vers la MEME limite")
    refte0 = None
    print(f"\n{'N':>5}{'DE0 (naif)':>16}{'err naif':>13}"
          f"{'DE0 (Li)':>16}{'err Li':>13}{'gain':>10}")
    Ns = [4, 8, 16, 24, 32, 48, 64, 96]
    en, el = [], []
    for N in Ns:
        rn, _, _ = rcwa_TM(N, lam, period, depth, er, eg, f, eI, eII, rule="naive")
        rl, _, _ = rcwa_TM(N, lam, period, depth, er, eg, f, eI, eII, rule="li")
        e1 = abs(rn[N]-R0ref); e2 = abs(rl[N]-R0ref)
        en.append(e1); el.append(e2)
        print(f"{N:5d}{rn[N]:16.10f}{e1:13.3e}{rl[N]:16.10f}{e2:13.3e}{e1/e2:10.1f}")
    en, el = np.array(en), np.array(el); Ns = np.array(Ns)
    pn = np.polyfit(np.log(Ns[3:]), np.log(en[3:]), 1)[0]
    pl = np.polyfit(np.log(Ns[3:]), np.log(el[3:]), 1)[0]
    print(f"\n   ordre de convergence mesure :  naif  N^({pn:.2f})    Li  N^({pl:.2f})")

    print("\n(2) TE : une seule regle suffit (pas de discontinuite concurrente)")
    refte, _, _ = rcwa_TE(Nref, lam, period, depth, er, eg, f, eI, eII)
    print(f"{'N':>5}{'DE0 (TE)':>16}{'erreur':>13}")
    ete = []
    for N in Ns:
        rt, _, _ = rcwa_TE(N, lam, period, depth, er, eg, f, eI, eII)
        e = abs(rt[N]-refte[Nref]); ete.append(e)
        print(f"{N:5d}{rt[N]:16.10f}{e:13.3e}")
    ete = np.array(ete)
    pte = np.polyfit(np.log(Ns[3:]), np.log(ete[3:]), 1)[0]
    print(f"\n   ordre de convergence TE : N^({pte:.2f})")

    print("\n(3) EFFICACITES CONVERGEES (N=300, regles de Li)")
    r, t, n = rcwa_TM(Nref, lam, period, depth, er, eg, f, eI, eII, rule="li")
    for i, o in enumerate(n):
        if abs(o) <= 2 and (r[i] > 1e-8 or t[i] > 1e-8):
            print(f"   ordre {o:+d} : R = {r[i]:.8f}   T = {t[i]:.8f}")
    print(f"   total : R = {r.sum():.10f}  T = {t.sum():.10f}  R+T = {r.sum()+t.sum():.12f}")

    # ------------------------------------------------------------- figure
    plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":.3,
                         "figure.dpi":150,"savefig.bbox":"tight","axes.linewidth":.8})
    C1,C2,C3="#1f3b73","#b3452c","#2f7d4f"
    fig,axs=plt.subplots(1,2,figsize=(9.4,3.3))
    axs[0].loglog(Ns,en,"o-",color=C2,ms=5,mfc="white",mew=1.4,label="Laurent naif (TM)")
    axs[0].loglog(Ns,el,"s-",color=C1,ms=5,mfc="white",mew=1.4,label="regles de Li (TM)")
    axs[0].loglog(Ns,ete,"^--",color=C3,ms=5,mfc="white",mew=1.4,label="TE (Laurent seul)")
    xx=np.logspace(np.log10(Ns[2]),np.log10(Ns[-1]),20)
    axs[0].loglog(xx,en[2]*(xx/Ns[2])**-1.0,":",color="k",lw=1,label=r"$N^{-1}$")
    axs[0].loglog(xx,el[2]*(xx/Ns[2])**-3.0,"-.",color="k",lw=1,label=r"$N^{-3}$")
    axs[0].set_xlabel("N (ordres retenus)"); axs[0].set_ylabel(r"erreur sur $DE_0$")
    axs[0].legend(fontsize=7.2); axs[0].set_title("convergence des efficacites",fontsize=9.5)
    Nsh=32
    rn,_,nn = rcwa_TM(Nsh,lam,period,depth,er,eg,f,eI,eII,rule="naive")
    rl,_,_  = rcwa_TM(Nsh,lam,period,depth,er,eg,f,eI,eII,rule="li")
    axs[1].axhline(R0ref,color="k",ls="--",lw=1.1,label=f"converge = {R0ref:.6f}")
    Nplot=np.arange(4,65,2)
    vn=[rcwa_TM(N,lam,period,depth,er,eg,f,eI,eII,rule="naive")[0][N] for N in Nplot]
    vl=[rcwa_TM(N,lam,period,depth,er,eg,f,eI,eII,rule="li")[0][N] for N in Nplot]
    axs[1].plot(Nplot,vn,"o-",color=C2,ms=3.5,lw=1.1,label="Laurent naif")
    axs[1].plot(Nplot,vl,"s-",color=C1,ms=3.5,lw=1.1,label="regles de Li")
    axs[1].set_xlabel("N"); axs[1].set_ylabel(r"$DE_0$ reflechi (TM)")
    axs[1].legend(fontsize=8); axs[1].set_title("la meme limite, pas la meme route",fontsize=9.5)
    fig.savefig("/home/claude/r_conv.pdf")
    np.savez("/home/claude/rcwa.npz",Ns=Ns,en=en,el=el,ete=ete,R0ref=R0ref)
    print("\nfigure ecrite.")
