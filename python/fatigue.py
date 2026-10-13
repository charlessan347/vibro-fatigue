# -*- coding: utf-8 -*-
"""
Post-traitement fatigue de la poutre encastrée sous excitation de base harmonique.

Entrée : reponse_1g.csv (extrait de l'ODB par extract_odb.py)
         -> S11(f) à la racine pour 1 g, amortissement 2 %

Étapes :
  1. Courbe S-N estimée (Basquin + facteurs de Marin) pour un acier S355
  2. Correction de contrainte moyenne (Goodman) : poids propre
  3. Cas 1 : excitation entretenue à la résonance (dwell) -> durée de vie
  4. Cas 2 : balayage sinus 10-1000 Hz (essai de qualification) -> Miner
  5. Sensibilité au niveau d'excitation et à l'amortissement
Unités : mm, t, s, MPa
"""
import os
import numpy as np
import matplotlib.pyplot as plt

DOSSIER = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(DOSSIER, 'reponse_1g.csv')
ZETA_REF = 0.02          # amortissement utilisé dans Abaqus
G = 9810.0               # mm/s²


# ---------------------------------------------------------------------------
# 1. Courbe S-N
# ---------------------------------------------------------------------------
class CourbeSN:
    """Courbe S-N estimée à partir de Su (méthode Shigley).

    - à 1e3 cycles : S = 0,9 Su (flexion)
    - à 1e6 cycles : Se = ka * kb * 0,5 Su (limite d'endurance corrigée)
    - entre les deux : droite de Basquin  S = a * N^b  en log-log
    """

    def __init__(self, Su, Sy, ka, kb, N1=1e3, Ne=1e6):
        self.Su, self.Sy = Su, Sy
        self.Se = ka * kb * 0.5 * Su
        self.S1 = 0.9 * Su
        self.N1, self.Ne = N1, Ne
        self.b = np.log10(self.S1 / self.Se) / np.log10(N1 / Ne)
        self.a = self.S1 / N1 ** self.b

    def N(self, S, limite=True):
        """Cycles à rupture pour une amplitude S.
        limite=True  : Miner original, S < Se -> vie infinie
        limite=False : Miner élémentaire, Basquin prolongé sous Se (conservatif)
        """
        S = np.asarray(S, dtype=float)
        with np.errstate(divide='ignore'):
            N = (S / self.a) ** (1.0 / self.b)
        if limite:
            N = np.where(S < self.Se, np.inf, N)
        return N

    def __str__(self):
        return ("Courbe S-N : Se = %.0f MPa, S(1e3) = %.0f MPa, "
                "Basquin a = %.0f MPa, b = %.4f" % (self.Se, self.S1, self.a, self.b))


def goodman(Sa, Sm, Su):
    """Amplitude alternée équivalente (contrainte moyenne nulle)."""
    return Sa / (1.0 - Sm / Su)


# ---------------------------------------------------------------------------
# 2. Réponse harmonique (issue d'Abaqus)
# ---------------------------------------------------------------------------
class ReponseHarmonique:
    def __init__(self, chemin):
        d = np.loadtxt(chemin, delimiter=',', skiprows=1)
        self.f, self.s = d[:, 0], d[:, 1]
        i = np.argmax(self.s)
        self.f_pic, self.s_pic = self.f[i], self.s[i]

    def S(self, f, g=1.0):
        """Amplitude S11 (MPa) à la fréquence f pour un niveau g (réponse linéaire).
        Interpolation en log pour suivre les pics de résonance."""
        return g * np.exp(np.interp(np.log(f), np.log(self.f), np.log(self.s)))


# ---------------------------------------------------------------------------
# 3. Cas de chargement
# ---------------------------------------------------------------------------
def dwell(S_amp, f, sn, limite=True):
    """Excitation entretenue à la fréquence f : (cycles, secondes) à rupture."""
    N = sn.N(S_amp, limite)
    return N, N / f


def balayage(rep, sn, g, R_oct_min=1.0, f0=10.0, f1=1000.0, limite=True, n=20000):
    """Dommage de Miner pour UN balayage sinus logarithmique f0 -> f1.

    Balayage log : f(t) = f0 * 2^(R t / 60)  =>  df/dt = f ln2 R / 60
    Cycles dans [f, f+df] : dn = f dt = 60 / (ln2 R) df  (constant par Hz)
    Dommage : D = somme dn / N(S(f))
    """
    f = np.geomspace(f0, f1, n)
    S = rep.S(f, g)
    dn_df = 60.0 / (np.log(2) * R_oct_min)
    dD_df = dn_df / sn.N(S, limite)
    integ = getattr(np, "trapezoid", None) or np.trapz
    D = integ(dD_df, f)
    return D, f, dD_df


# ---------------------------------------------------------------------------
# 4. Script principal
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # --- Matériau : acier S355, surface usinée -------------------------------
    Su, Sy = 510.0, 355.0
    ka = 4.51 * Su ** -0.265                  # état de surface usiné (Shigley)
    de = 0.808 * np.sqrt(20 * 5)              # diamètre équivalent section rect.
    kb = 1.24 * de ** -0.107                  # effet d'échelle
    sn = CourbeSN(Su, Sy, ka, kb)
    print(sn)
    print("ka = %.3f, kb = %.3f" % (ka, kb))

    rep = ReponseHarmonique(CSV)
    print("\nPic Abaqus (1 g, zeta = 2 %%) : %.1f MPa à %.2f Hz" % (rep.s_pic, rep.f_pic))

    # --- Contrainte moyenne : poids propre ----------------------------------
    rho, A, L, b_, h_ = 7.85e-9, 100.0, 300.0, 20.0, 5.0
    M_racine = rho * A * G * L ** 2 / 2       # moment d'encastrement (N.mm)
    Sm = M_racine / (b_ * h_ ** 2 / 6)
    print("Contrainte moyenne (poids propre) : %.1f MPa -> Goodman : x%.3f"
          % (Sm, 1.0 / (1.0 - Sm / Su)))

    # --- Niveau limite pour une vie infinie ---------------------------------
    g_lim = sn.Se / goodman(rep.s_pic, Sm, Su)
    g_plast = sn.Sy / rep.s_pic
    print("\nNiveau limite vie infinie (dwell résonance) : %.2f g" % g_lim)
    print("Niveau de début de plastification         : %.2f g" % g_plast)

    # --- Cas 1 : dwell à la résonance ---------------------------------------
    print("\n--- Dwell à la résonance (%.1f Hz) ---" % rep.f_pic)
    print("  g    S_a (MPa)   N (cycles)     durée")
    for g in [1, 2, 2.5, 3, 3.5]:
        Sa = goodman(g * rep.s_pic, Sm, Su)
        N, t = dwell(Sa, rep.f_pic, sn)
        duree = "infinie" if np.isinf(t) else "%.1f min" % (t / 60)
        flag = "  (> Sy : hors linéaire)" if g * rep.s_pic > Sy else ""
        print("  %.1f  %7.1f   %12.3g   %s%s" % (g, Sa, N, duree, flag))

    # --- Cas 2 : balayage sinus ---------------------------------------------
    g_bal, R = 3.0, 1.0
    D, f, dD = balayage(rep, sn, g_bal, R)
    D_el, _, _ = balayage(rep, sn, g_bal, R, limite=False)
    duree_bal = np.log2(1000 / 10) / R
    print("\n--- Balayage 10-1000 Hz à %.1f g, %.0f oct/min (%.1f min/balayage) ---"
          % (g_bal, R, duree_bal))
    print("  Dommage par balayage (Miner original)   : %.3e -> %.0f balayages"
          % (D, 1 / D if D > 0 else np.inf))
    print("  Dommage par balayage (Miner élémentaire): %.3e -> %.0f balayages"
          % (D_el, 1 / D_el))

    # --- Sensibilité à l'amortissement (dwell, 3 g) -------------------------
    print("\n--- Sensibilité à l'amortissement (dwell résonance, 3 g) ---")
    for z in [0.01, 0.02, 0.03, 0.05]:
        S_brut = 3 * rep.s_pic * ZETA_REF / z      # pic de résonance en 1/zeta
        Sa = goodman(S_brut, Sm, Su)
        N, t = dwell(Sa, rep.f_pic, sn)
        if S_brut > Sy:
            duree = "hors domaine (S > Sy, calcul linéaire non valide)"
        else:
            duree = "infinie" if np.isinf(t) else "%.1f min" % (t / 60)
        print("  zeta = %.0f %%  S_a = %6.1f MPa  -> %s" % (100 * z, Sa, duree))

    # --- Figures -------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

    ax[0].loglog(rep.f, rep.s, 'o-', ms=3)
    ax[0].set(xlabel='Fréquence (Hz)', ylabel='S11 racine (MPa)',
              title='Fonction de transfert en contrainte (1 g)')
    ax[0].grid(True, which='both', alpha=0.3)

    Nn = np.logspace(3, 8, 200)
    Sn = np.where(Nn < sn.Ne, sn.a * Nn ** sn.b, sn.Se)
    ax[1].loglog(Nn, Sn, 'k-', label='Courbe S-N estimée')
    for g in [2.5, 3, 3.5]:
        Sa = goodman(g * rep.s_pic, Sm, Su)
        ax[1].plot(sn.N(Sa), Sa, 'o', label='Dwell %.1f g' % g)
    ax[1].axhline(Sy, ls='--', c='grey', lw=0.8)
    ax[1].text(1.2e3, Sy * 1.03, 'Sy', color='grey')
    ax[1].set(xlabel='Cycles à rupture N', ylabel='Amplitude (MPa)', title='Courbe S-N')
    ax[1].legend(); ax[1].grid(True, which='both', alpha=0.3)

    _, f_el, dD_el = balayage(rep, sn, g_bal, R, limite=False)
    ax[2].semilogx(f_el, dD_el, label='Miner élémentaire')
    ax[2].semilogx(f, dD, label='Miner original')
    ax[2].set(xlabel='Fréquence (Hz)', ylabel='dD/df (1/Hz)',
              title='Densité de dommage, balayage %.0f g' % g_bal, xlim=(30, 70))
    ax[2].legend(); ax[2].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(DOSSIER, 'fatigue.png'), dpi=150)
    plt.show()
