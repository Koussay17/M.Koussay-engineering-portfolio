#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Instabilite double faisceau : theorie lineaire de Vlasov-Poisson
et verification par un code particulaire 1D-1V ecrit independamment.
Unites : omega_p = m = e = eps0 = n0 = 1.
"""
import numpy as np
from scipy.special import wofz
from scipy.optimize import brentq
from scipy.integrate import quad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260727)

# ====================================================================
#  1. THEORIE LINEAIRE
# ====================================================================
def cold_gamma(a):
    """a = k v0 ; taux de croissance des faisceaux froids."""
    Xm = (2*a**2 + 1 - np.sqrt(1 + 8*a**2))/2
    return np.sqrt(-Xm) if Xm < 0 else 0.0

def Zfun(z):
    return 1j*np.sqrt(np.pi)*wofz(z)

def eps_im(k, g, v0, vt):
    """eps(k, i gamma) ; reel par symetrie."""
    s = 0.0
    for al, vj in ((0.5, v0), (0.5, -v0)):
        zj = (1j*g - k*vj)/(np.sqrt(2)*k*vt)
        s += al*(1.0 + zj*Zfun(zj))
    return (1.0 + s/(k**2*vt**2)).real

def warm_gamma(k, v0, vt, gmax=1.0):
    """Racine de eps(k,i gamma)=0 sur gamma>0, ou 0 si stable."""
    gs = np.linspace(1e-6, gmax, 4000)
    vals = np.array([eps_im(k, g, v0, vt) for g in gs])
    sign = np.sign(vals)
    idx = np.where(sign[:-1] != sign[1:])[0]
    if len(idx) == 0:
        return 0.0
    i = idx[-1]
    return brentq(lambda g: eps_im(k, g, v0, vt), gs[i], gs[i+1], xtol=1e-13)

def penrose(v0, vt):
    f0 = lambda v: (np.exp(-(v-v0)**2/(2*vt**2))
                    + np.exp(-(v+v0)**2/(2*vt**2)))/(2*np.sqrt(2*np.pi)*vt)
    f0m = f0(0.0); L = v0 + 14*vt
    g = lambda v: (f0(v) - f0m)/v**2
    return quad(g, -L, 0, limit=400)[0] + quad(g, 0, L, limit=400)[0]

# ====================================================================
#  2. CODE PARTICULAIRE (PIC) 1D-1V ELECTROSTATIQUE
# ====================================================================
class PIC1D:
    """CIC + Poisson spectral + saute-mouton. Ions fixes neutralisants."""
    def __init__(self, L, Ng, Nx, Nv, v0, vt, dt, amp=1e-6, nmodes=12, seed=7):
        """Chargement 'quiet start' tensoriel : Nx positions x Nv vitesses par
        faisceau, donc densite initiale exactement uniforme (bruit de grenaille
        nul). Une perturbation controlee d'amplitude amp ensemence tous les
        modes m = 1..nmodes."""
        N = 2*Nx*Nv
        self.L, self.Ng, self.N, self.dt = L, Ng, N, dt
        self.dx = L/Ng
        self.w  = L/N
        from scipy.special import ndtri
        xa = L*(np.arange(Nx)+0.5)/Nx
        vb = ndtri((np.arange(Nv)+0.5)/Nv)*vt
        X  = np.tile(xa, Nv)
        V  = np.repeat(vb, Nx)
        self.x = np.concatenate([X, X])
        self.v = np.concatenate([V+v0, V-v0])
        r = np.random.default_rng(seed)
        ph = r.random(nmodes)*2*np.pi
        dxp = np.zeros_like(self.x)
        for m in range(1, nmodes+1):
            km = 2*np.pi*m/L
            dxp += amp*L*np.sin(km*self.x + ph[m-1])
        self.x = (self.x + dxp) % L
        kk = 2*np.pi*np.fft.rfftfreq(Ng, d=self.dx)
        self.k = kk
        self.k2 = np.where(kk == 0, 1.0, kk**2)

    def deposit(self):
        xg = self.x/self.dx
        i0 = np.floor(xg).astype(np.int64) % self.Ng
        f  = xg - np.floor(xg)
        rho = np.zeros(self.Ng)
        np.add.at(rho, i0, (1-f))
        np.add.at(rho, (i0+1) % self.Ng, f)
        rho *= -self.w/self.dx          # densite de charge electronique
        return rho + 1.0                # + fond ionique fixe (n0 = 1)

    def field(self, rho):
        rk = np.fft.rfft(rho)
        Ek = -1j*rk/np.where(self.k == 0, 1.0, self.k)
        Ek[0] = 0.0
        return np.fft.irfft(Ek, n=self.Ng), rk

    def gather(self, E):
        xg = self.x/self.dx
        i0 = np.floor(xg).astype(np.int64) % self.Ng
        f  = xg - np.floor(xg)
        return E[i0]*(1-f) + E[(i0+1) % self.Ng]*f

    def run(self, tmax, snap_times=()):
        nst = int(round(tmax/self.dt))
        rho = self.deposit(); E, _ = self.field(rho)
        a = -self.gather(E)
        self.v -= 0.5*self.dt*a                      # demi-pas arriere
        t = 0.0
        hist = dict(t=[], Ek=[], WE=[], WK=[], snaps=[])
        nk = len(self.k)
        for n in range(nst):
            rho = self.deposit(); E, rk = self.field(rho)
            a = -self.gather(E)
            self.v += self.dt*a
            self.x = (self.x + self.dt*self.v) % self.L
            t += self.dt
            amp = np.abs(rk)/self.Ng               # amplitude des modes de rho
            hist['t'].append(t); hist['Ek'].append(amp.copy())
            hist['WE'].append(0.5*np.sum(E**2)*self.dx)
            hist['WK'].append(0.5*self.w*np.sum(self.v**2))
            for ts in snap_times:
                if abs(t-ts) < 0.5*self.dt:
                    sel = rng.choice(self.N, size=min(60000, self.N), replace=False)
                    hist['snaps'].append((t, self.x[sel].copy(), self.v[sel].copy()))
        for key in ('t','WE','WK'):
            hist[key] = np.array(hist[key])
        hist['Ek'] = np.array(hist['Ek'])
        return hist

# ====================================================================
#  3. EXECUTION
# ====================================================================
if __name__ == "__main__":
    v0, vt = 1.0, 0.15
    L, Ng, dt, tmax = 20*np.pi, 256, 0.05, 75.0
    Nx, Nv = 2048, 256
    N = 2*Nx*Nv
    kmodes = 2*np.pi*np.arange(1, 13)/L

    print("="*70); print("THEORIE LINEAIRE"); print("="*70)
    aa = np.linspace(1e-5, 0.99999, 400001)
    gg = np.array([cold_gamma(x) for x in aa]); i = np.argmax(gg)
    print(f"  froid  : a_max = {aa[i]:.6f} (exact sqrt(3/8) = {np.sqrt(3/8):.6f})")
    print(f"           gamma_max = {gg[i]:.6f} (exact 1/(2sqrt2) = {1/(2*np.sqrt(2)):.6f})")
    rc = brentq(lambda r: penrose(r, 1.0), 0.6, 2.0, xtol=1e-12)
    print(f"  Penrose: seuil v0/vt = {rc:.8f}")
    print(f"  chaud  : v0 = {v0}, vt = {vt}  (v0/vt = {v0/vt:.2f})")
    kfine = np.linspace(0.02, 1.3, 300)
    gw = np.array([warm_gamma(k, v0, vt) for k in kfine])
    j = np.argmax(gw)
    print(f"           k_max = {kfine[j]:.5f}, gamma_max = {gw[j]:.6f}")
    gth = np.array([warm_gamma(k, v0, vt) for k in kmodes])
    for k, g in zip(kmodes, gth):
        print(f"           k = {k:.4f} -> gamma_th = {g:.6f}")

    print(); print("="*70)
    print(f"PIC : L={L:.4f}, Ng={Ng}, N={N:,}, dt={dt}, tmax={tmax}")
    print("="*70)
    sim = PIC1D(L, Ng, Nx, Nv, v0, vt, dt)
    h = sim.run(tmax, snap_times=(0.0, 26.0, 34.0, 55.0))
    t, Ek, WE, WK = h['t'], h['Ek'], h['WE'], h['WK']

    # ajustement du taux de croissance sur la phase lineaire
    print(f"{'k':>8}{'gamma_PIC':>12}{'gamma_th':>12}{'ecart %':>10}   fenetre")
    gpic = []
    for m, (k, gt) in enumerate(zip(kmodes, gth), start=1):
        y = np.log(Ek[:, m] + 1e-300)
        if gt > 0.02:
            t1, t2 = 6.0, min(6.0 + 6.0/gt, 24.0)
        else:
            t1, t2 = 6.0, 24.0
        w = (t > t1) & (t < t2)
        p = np.polyfit(t[w], y[w], 1)
        gpic.append(p[0])
        err = 100*(p[0]-gt)/gt if gt > 1e-3 else np.nan
        print(f"{k:8.4f}{p[0]:12.5f}{gt:12.5f}{err:10.2f}   [{t1:.0f},{t2:.0f}]")
    gpic = np.array(gpic)

    # saturation et piegeage
    isat = np.argmax(WE)
    Emax = np.sqrt(2*WE[isat]/L)*np.sqrt(2)      # amplitude ~ sqrt(2<E^2>)
    kdom = kmodes[np.argmax(gth)]
    Etrap = gth.max()**2/kdom
    print(f"\n  saturation : t_sat = {t[isat]:.2f}, W_E^max = {WE[isat]:.5f}")
    print(f"  amplitude E_0 mesuree ~ {Emax:.4f}")
    print(f"  estimation par piegeage gamma^2/k = {Etrap:.4f}  (rapport {Emax/Etrap:.3f})")
    Etot = WE + WK
    print(f"  conservation de l'energie : dE/E0 = {(Etot.max()-Etot.min())/Etot[0]:.3e}")

    np.savez("/home/claude/pic.npz", t=t, Ek=Ek, WE=WE, WK=WK, kmodes=kmodes,
             gth=gth, gpic=gpic, kfine=kfine, gw=gw, L=L)
    import pickle
    with open("/home/claude/snaps.pkl","wb") as f: pickle.dump(h['snaps'], f)
    print("\ndonnees sauvegardees.")
