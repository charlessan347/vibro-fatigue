import numpy as np

# Unités : mm, t, s, MPa
E, nu, rho = 210000, 0.3, 7.85e-9
L, b, h = 300, 20, 5
G = E / (2*(1 + nu))
A = b * h
M = rho * A * L
I_faible = b * h**3 / 12
I_forte  = h * b**3 / 12

# --- Flexion (Euler-Bernoulli) ---
beta_L = np.array([1.875, 4.694, 7.855, 10.996, 14.137])
def freq_flexion(I):
    return beta_L**2 / (2*np.pi) * np.sqrt(E*I / (rho*A*L**4))

print("Flexion faible (Hz):", freq_flexion(I_faible).round(1))
print("Flexion forte  (Hz):", freq_flexion(I_forte).round(1))

# --- Torsion (Saint-Venant) ---
J  = 0.281 * b * h**3                 # coeff. tabulé pour b/h = 4
Ip = b * h * (b**2 + h**2) / 12
n  = np.arange(1, 3)
f_tor = (2*n - 1) / (4*L) * np.sqrt(G*J / (rho*Ip))
print("Torsion        (Hz):", f_tor.round(1))

# --- Réponse harmonique au mode 1 (excitation de base) ---
a0, zeta = 9810, 0.02                 # 1 g en mm/s², amortissement 2 %
w1 = 2*np.pi * freq_flexion(I_faible)[0]
sigma1 = 0.7341                       # coeff. de forme du mode 1
Gamma = np.sqrt(M) * 2*sigma1 / beta_L[0]   # facteur de participation
q = Gamma * a0 / (w1**2 * 2*zeta)           # amplitude modale à la résonance
w_bout = 2 / np.sqrt(M) * q                 # déplacement relatif en bout
kappa  = w_bout * (beta_L[0] / L)**2        # courbure à l'encastrement
sig    = E * kappa * h / 2                  # contrainte de flexion en peau

print(f"Gamma mode 1      : {Gamma:.4f}   (Abaqus : 0.0120)")
print(f"Déplacement bout  : {w_bout:.2f} mm")
print(f"Contrainte racine : {sig:.1f} MPa")