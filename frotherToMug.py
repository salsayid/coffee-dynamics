import numpy as np
import matplotlib.pyplot as plt

# parameters
R1 = 0.0105       # frother disc radius - m
R2 = 0.035        # mug inner radius - m
Omega1 = 1257     # angular velocity - rad/s

# coefficients from boundary conditions
A = (-Omega1 * R1**2) / (R2**2 - R1**2)
B = (Omega1 * R1**2 * R2**2) / (R2**2 - R1**2)

# radial grid
r = np.linspace(R1, R2, 500)

# velocity profile
u = A * r + B / r

# plotting
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(r * 1000, u, color='#2c5f8a', linewidth=2)

ax.set_xlabel('Radial position  r  [mm]', fontsize=12)
ax.set_ylabel('Azimuthal velocity  $u_\\theta$  [m/s]', fontsize=12)
ax.set_title('Laminar Velocity Profile — Steady State', fontsize=13)

ax.axvline(R1 * 1000, color='gray', linewidth=0.8, linestyle='--', label='Frother disc $R_1$')
ax.axvline(R2 * 1000, color='gray', linewidth=0.8, linestyle=':', label='Mug wall $R_2$')

ax.annotate('Frother disc\n$u_\\theta = \\Omega_1 R_1 = 13.2$ m/s',
            xy=(R1 * 1000, Omega1 * R1),
            xytext=(14, 10),
            arrowprops=dict(arrowstyle='->', color='black'),
            fontsize=9)

ax.annotate('Mug wall\n$u_\\theta = 0$',
            xy=(R2 * 1000, 0),
            xytext=(28, 2),
            arrowprops=dict(arrowstyle='->', color='black'),
            fontsize=9)

ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('laminar_profile.png', dpi=200, bbox_inches='tight')
plt.show()

print(f"A = {A:.2f} s⁻¹")
print(f"B = {B:.6f} m²/s")
print(f"u_theta(R1) = {A*R1 + B/R1:.2f} m/s  (should be {Omega1*R1:.2f})")
print(f"u_theta(R2) = {A*R2 + B/R2:.4f} m/s  (should be 0.0000)")
