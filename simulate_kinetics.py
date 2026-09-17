import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. PARAMETERS & INITIAL CONCENTRATIONS (nM)
# ------------------------------------------------------------------------------
# Rates (typical CHA bi-molecular & unimolecular rate constants)
k1 = 1e5    # M^-1 s^-1 (Target + H1 binding)
k1_r = 0.1  # s^-1 (Unbinding)
k2 = 1e5    # M^-1 s^-1 (H1:Target + H2 displacement)
k_leak = 10 # M^-1 s^-1 (Spontaneous H1 + H2 leakage)

# Initial concentrations in Molar (from config)
H1_0 = 1e-7
H2_0 = 5e-8
Target_0 = 2.5e-8

# ------------------------------------------------------------------------------
# 2. ODE SYSTEM DEFINITION
# ------------------------------------------------------------------------------
# States: y = [Target, H1, H2, H1:Target, Complex_3Strand, Leak_H1H2]
def cha_system(t, y):
    # Enforce non-negative bounds to prevent solver overshoots
    T, H1, H2, H1T, C3, Leak = [max(0.0, val) for val in y]
    
    dT = -k1 * T * H1 + k1_r * H1T + k2 * H1T * H2
    dH1 = -k1 * T * H1 + k1_r * H1T - k_leak * H1 * H2
    dH2 = -k2 * H1T * H2 - k_leak * H1 * H2
    dH1T = k1 * T * H1 - k1_r * H1T - k2 * H1T * H2
    dC3 = k2 * H1T * H2
    dLeak = k_leak * H1 * H2
    
    return [dT, dH1, dH2, dH1T, dC3, dLeak]

# Time span: 0 to 7200 seconds (2 hours)
t_span = (0, 7200)
t_eval = np.linspace(0, 7200, 500)
y0 = [Target_0, H1_0, H2_0, 0, 0, 0]

# ------------------------------------------------------------------------------
# 3. RUN SIMULATION & PLOT
# ------------------------------------------------------------------------------
# Using Radau method designed for stiff chemical reaction systems
sol = solve_ivp(cha_system, t_span, y0, t_eval=t_eval, method='Radau')

signal = sol.y[4] * 1e9  # Catalytic complex (nM)
leak = sol.y[5] * 1e9    # Background leak (nM)
time_min = sol.t / 60

snr = signal[-1] / (leak[-1] + 1e-12)

print(f"=== KINETIC SIMULATION RESULTS (2 Hours) ===")
print(f"Final Signal (Active Complex) : {signal[-1]:.2f} nM")
print(f"Final Background (Leak)       : {leak[-1]:.2f} nM")
print(f"Signal-to-Noise Ratio (SNR)   : {snr:.2f}x")

plt.figure(figsize=(8, 5))
plt.plot(time_min, signal, label="Target-Triggered Signal (C3)", color="green", linewidth=2)
plt.plot(time_min, leak, label="Un-triggered Background Leak", color="red", linestyle="--", linewidth=2)
plt.xlabel("Time (minutes)")
plt.ylabel("Concentration (nM)")
plt.title("CHA Reaction Kinetics (Candidate 28)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("cha_kinetics.png")
print("\nPlot saved as 'cha_kinetics.png'.")