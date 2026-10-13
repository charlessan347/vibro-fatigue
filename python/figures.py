# -*- coding: utf-8 -*-
"""
Courbe de réponse en fréquence (1 g, zeta = 2 %) à partir de reponse_1g.csv.
Sortie : images/reponse_frequence.png
Lancement : dans VSCode, bouton ▶ (après extract_odb.py)
"""
import os
import numpy as np
import matplotlib.pyplot as plt

ICI = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(ICI, 'reponse_1g.csv')
SORTIE = os.path.join(ICI, '..', 'images', 'reponse_frequence.png')
os.makedirs(os.path.dirname(SORTIE), exist_ok=True)

f, s11, u3 = np.loadtxt(CSV, delimiter=',', skiprows=1, unpack=True)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ax1.loglog(f, u3, color='#1d5f8a')
ax1.set_ylabel('U3 en bout (mm)')
ax1.set_title('Réponse à une excitation de base de 1 g (amortissement 2 %)')

ax2.loglog(f, s11, color='#b5502a')
ax2.set_ylabel('S11 à l\'encastrement (MPa)')
ax2.set_xlabel('Fréquence (Hz)')

for ax in (ax1, ax2):
    ax.grid(True, which='both', alpha=0.3)
    for fr, nom in [(46.6, 'mode 1'), (291.8, 'mode 3'), (815.9, 'mode 4')]:
        ax.axvline(fr, color='grey', lw=0.8, ls='--')
ymax = ax1.get_ylim()[1]
for fr, nom in [(46.6, 'flexion 1'), (291.8, 'flexion 2'), (815.9, 'flexion 3')]:
    ax1.text(fr * 1.05, ymax * 0.3, nom, color='grey', fontsize=9)

fig.tight_layout()
fig.savefig(SORTIE, dpi=150)
print('Figure écrite :', os.path.abspath(SORTIE))
plt.show()