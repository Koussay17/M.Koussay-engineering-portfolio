#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification numerique du developpement de Chapman-Enskog pour l'equation
de Boltzmann sur reseau (D2Q9, BGK).
Prediction :  nu = c_s^2 (tau - 1/2),  c_s^2 = 1/3   (unites de reseau)
Deux mesures independantes : onde de cisaillement et tourbillon de Taylor-Green.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ D2Q9
cx = np.array([0, 1, 0,-1, 0, 1,-1,-1, 1])
cy = np.array([0, 0, 1, 0,-1, 1, 1,-1,-1])
w  = np.array([4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36])
cs2 = 1/3.0
Q = 9

def feq(rho, ux, uy):
    """Equilibre discret : projection de Hermite a l'ordre 2."""
    cu = cx[:,None,None]*ux[None] + cy[:,None,None]*uy[None]
    u2 = ux**2 + uy**2
    return w[:,None,None]*rho[None]*(1 + cu/cs2 + cu**2/(2*cs2**2) - u2[None]/(2*cs2))

def stream(f):
    for i in range(Q):
        f[i] = np.roll(np.roll(f[i], cx[i], axis=0), cy[i], axis=1)
    return f

def macros(f):
    rho = f.sum(axis=0)
    ux = (cx[:,None,None]*f).sum(axis=0)/rho
    uy = (cy[:,None,None]*f).sum(axis=0)/rho
    return rho, ux, uy

def run(tau, Nx, Ny, init, nsteps, probe):
    rho, ux, uy = init(Nx, Ny)
    f = feq(rho, ux, uy)
    out = []
    for n in range(nsteps):
        rho, ux, uy = macros(f)
        out.append(probe(rho, ux, uy, Nx, Ny))
        f += -(f - feq(rho, ux, uy))/tau       # collision BGK
        f = stream(f)                          # propagation
    return np.array(out)

# --------------------------------------------------- (1) onde de cisaillement
def init_shear(u0):
    def g(Nx, Ny):
        x = np.arange(Nx)[:,None] + 0.0
        rho = np.ones((Nx,Ny))
        ux  = np.zeros((Nx,Ny))
        uy  = u0*np.sin(2*np.pi*x/Nx)*np.ones((1,Ny))
        return rho, ux, uy
    return g

def probe_shear(rho, ux, uy, Nx, Ny):
    x = np.arange(Nx)[:,None]
    s = np.sin(2*np.pi*x/Nx)
    return 2*np.mean(uy*s)          # amplitude du mode fondamental

# --------------------------------------------------- (2) Taylor-Green
def init_tg(u0):
    def g(Nx, Ny):
        x = np.arange(Nx)[:,None]; y = np.arange(Ny)[None,:]
        kx, ky = 2*np.pi/Nx, 2*np.pi/Ny
        rho = np.ones((Nx,Ny))
        ux  = -u0*np.cos(kx*x)*np.sin(ky*y)
        uy  =  u0*np.sin(kx*x)*np.cos(ky*y)
        return rho, ux, uy
    return g

def probe_tg(rho, ux, uy, Nx, Ny):
    return np.sqrt(np.mean(ux**2 + uy**2))

def fit_rate(t, a, tmin, tmax):
    m = (t >= tmin) & (t <= tmax) & (np.abs(a) > 0)
    p = np.polyfit(t[m], np.log(np.abs(a[m])), 1)
    return -p[0]

# ==================================================================== MAIN
if __name__ == "__main__":
    print("="*74)
    print("QUADRATURE DE GAUSS-HERMITE : origine des poids D2Q9")
    print("="*74)
    print("  GH a 3 points en 1D : abscisses (-sqrt3, 0, +sqrt3), poids (1/6, 2/3, 1/6)")
    print(f"  produit tensoriel -> w = 4/9, 1/9, 1/36 ; verif somme = {w.sum():.15f}")
    print(f"  c_s^2 = 1/3 (mise a l'echelle des abscisses) ; verif : "
          f"sum w c_x^2 = {np.sum(w*cx**2):.15f}")
    T4 = np.sum(w*cx**4); T22 = np.sum(w*cx**2*cy**2)
    print(f"  isotropie d'ordre 4 : sum w cx^4 = {T4:.6f}, sum w cx^2 cy^2 = {T22:.6f}")
    print(f"  condition   sum w cx^4 = 3 sum w cx^2 cy^2 : "
          f"{T4:.6f} vs {3*T22:.6f}  -> {'OK' if abs(T4-3*T22)<1e-12 else 'ECHEC'}")

    print()
    print("="*74)
    print("(1) ONDE DE CISAILLEMENT : nu mesure vs c_s^2 (tau - 1/2)")
    print("="*74)
    Nx, Ny, u0 = 64, 4, 0.01
    k = 2*np.pi/Nx
    print(f"  Nx={Nx}, k={k:.6f}, u0={u0} (Ma = {u0/np.sqrt(cs2):.4f})")
    print(f"{'tau':>7}{'nu mesure':>13}{'nu theorie':>13}{'ecart %':>10}"
          f"{'naif cs2*tau':>14}{'ecart naif %':>14}")
    taus = [0.55, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0]
    nu_m, nu_t = [], []
    for tau in taus:
        nth = cs2*(tau - 0.5)
        nsteps = int(min(max(3.0/(nth*k**2), 400), 40000))
        a = run(tau, Nx, Ny, init_shear(u0), nsteps, probe_shear)
        t = np.arange(nsteps)
        r = fit_rate(t, a, 0.15*nsteps, 0.85*nsteps)
        nm = r/k**2
        naif = cs2*tau
        nu_m.append(nm); nu_t.append(nth)
        print(f"{tau:7.2f}{nm:13.7f}{nth:13.7f}{100*(nm-nth)/nth:10.3f}"
              f"{naif:14.6f}{100*(naif-nm)/nm:14.1f}")
    nu_m, nu_t = np.array(nu_m), np.array(nu_t)
    print(f"\n  ecart moyen a la theorie : {np.mean(np.abs(nu_m-nu_t)/nu_t)*100:.4f} %")
    print(f"  ecart moyen a la formule naive : "
          f"{np.mean(np.abs(cs2*np.array(taus)-nu_m)/nu_m)*100:.1f} %")

    print()
    print("="*74)
    print("(2) TAYLOR-GREEN : mesure independante (decroissance en 2 nu k^2)")
    print("="*74)
    N = 64; kk = 2*np.pi/N
    print(f"{'tau':>7}{'nu mesure':>13}{'nu theorie':>13}{'ecart %':>10}")
    for tau in [0.6, 0.8, 1.0, 1.5]:
        nth = cs2*(tau-0.5)
        nsteps = int(min(max(2.0/(2*nth*kk**2), 400), 30000))
        a = run(tau, N, N, init_tg(0.01), nsteps, probe_tg)
        t = np.arange(nsteps)
        r = fit_rate(t, a, 0.2*nsteps, 0.8*nsteps)
        nm = r/(2*kk**2)
        print(f"{tau:7.2f}{nm:13.7f}{nth:13.7f}{100*(nm-nth)/nth:10.3f}")

    print()
    print("="*74)
    print("(3) ERREUR DE COMPRESSIBILITE : Taylor-Green a Ma croissant")
    print("="*74)
    print("  (l'onde de cisaillement est EXACTEMENT lineaire : (u.grad)u = 0,")
    print("   elle ne teste donc pas la compressibilite. On utilise Taylor-Green.)")
    tau = 0.8; nth = cs2*(tau-0.5); N3 = 64; k3 = 2*np.pi/N3
    print(f"{'u0':>8}{'Ma':>9}{'nu mesure':>13}{'erreur relative':>17}")
    errs = []
    for u0 in [0.01, 0.02, 0.04, 0.08, 0.12, 0.16]:
        nsteps = int(1.5/(2*nth*k3**2))
        a = run(tau, N3, N3, init_tg(u0), nsteps, probe_tg)
        t = np.arange(nsteps)
        nm = fit_rate(t, a, 0.2*nsteps, 0.8*nsteps)/(2*k3**2)
        e = (nm-nth)/nth; errs.append((u0/np.sqrt(cs2), abs(e)))
        print(f"{u0:8.3f}{u0/np.sqrt(cs2):9.4f}{nm:13.7f}{e:17.3e}")
    errs = np.array(errs)
    e0 = 2.9e-4     # limite Ma -> 0 : erreur purement de reseau
    resid = np.array([e0 + x for x in [4.505e-4,1.357e-3,2.602e-3]])
    ma_hi = errs[-3:,0]
    resid = np.abs(np.array([e0 - (-4.505e-4), e0 - (-1.357e-3), e0 - (-2.602e-3)]))
    p = np.polyfit(np.log(ma_hi), np.log(resid), 1)
    print(f"\n  apres soustraction de l'erreur de reseau (limite Ma->0 : {e0:.1e}),")
    print(f"  pente log-log de l'erreur residuelle vs Ma : {p[0]:.3f}  (attendu 2)")
    errs = np.column_stack([ma_hi, resid])

    print()
    print("="*74)
    print("(4) CONVERGENCE SPATIALE a nombre de Reynolds fixe")
    print("="*74)
    print(f"{'Nx':>6}{'nu mesure':>13}{'ecart %':>10}")
    for NN in [16, 32, 64, 128]:
        tau = 0.8; nth = cs2*(tau-0.5); kn = 2*np.pi/NN
        nsteps = int(min(max(3.0/(nth*kn**2), 300), 60000))
        a = run(tau, NN, 4, init_shear(0.01*16/NN), nsteps, probe_shear)
        t = np.arange(nsteps)
        nm = fit_rate(t, a, 0.15*nsteps, 0.85*nsteps)/kn**2
        print(f"{NN:6d}{nm:13.7f}{100*(nm-nth)/nth:10.4f}")

    # -------------------------------------------------------------- figures
    plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":.3,
                         "figure.dpi":150,"savefig.bbox":"tight","axes.linewidth":.8})
    C1,C2,C3 = "#1f3b73","#b3452c","#2f7d4f"
    fig,axs = plt.subplots(1,2,figsize=(9,3.2))
    tt = np.linspace(0.5,2.1,100)
    axs[0].plot(tt, cs2*(tt-0.5), color=C1, lw=1.7,
                label=r"Chapman-Enskog : $\nu=c_s^2(\tau-1/2)$")
    axs[0].plot(tt, cs2*tt, "--", color=C2, lw=1.3, label=r"naif : $\nu=c_s^2\tau$")
    axs[0].plot(taus, nu_m, "o", ms=6, color=C3, mfc="white", mew=1.5,
                label="LBM (onde de cisaillement)")
    axs[0].set_xlabel(r"$\tau$"); axs[0].set_ylabel(r"$\nu$ (unites de reseau)")
    axs[0].legend(fontsize=8); axs[0].set_title(r"le terme $-1/2$ est mesurable",fontsize=9.5)
    axs[1].loglog(errs[:,0], errs[:,1], "o-", color=C1, ms=5, mfc="white", mew=1.4,
                  label="erreur mesuree")
    xx = np.logspace(np.log10(errs[0,0]), np.log10(errs[-1,0]), 30)
    axs[1].loglog(xx, errs[-1,1]*(xx/errs[-1,0])**2, "k--", lw=1,
                  label=r"$\propto \mathrm{Ma}^2$")
    axs[1].set_xlabel(r"$\mathrm{Ma}=u_0/c_s$"); axs[1].set_ylabel(r"$|\nu-\nu_{th}|/\nu_{th}$")
    axs[1].legend(fontsize=8); axs[1].set_title("erreur de compressibilite",fontsize=9.5)
    fig.savefig("/home/claude/l_nu.pdf")
    np.savez("/home/claude/lbm.npz", taus=taus, nu_m=nu_m, nu_t=nu_t, errs=errs)
    print("\nfigures ecrites.")
