"""Single source of truth for every geometric and physical parameter (SI units).

Values marked ESTIMATED are not measured; change them here and rebuild.
"""
import math

# --- fluid: water at 85 C -------------------------------------------------
RHO = 968.6            # kg/m^3
NU = 3.3e-7            # m^2/s
MU = RHO * NU          # Pa s
SIGMA = 0.0618         # N/m, water/air surface tension
RHO_AIR = 1.0          # kg/m^3, hot humid air
MU_AIR = 2.1e-5        # Pa s
G = 9.81               # m/s^2

# --- frother ----------------------------------------------------------------
RPM = 12000.0
OMEGA = RPM * 2 * math.pi / 60      # 1256.6 rad/s
R1 = 0.0105                         # disc / coil outer radius (README)

# coil whisk head, all ESTIMATED
COIL_RING_R = 0.0090        # torus centerline radius
COIL_D = 0.0030             # helical coil outer diameter (R1 = ring_r + coil_d/2)
COIL_WIRE_D = 0.0004        # spring wire diameter
COIL_PITCH = 0.0008         # pitch along the ring centerline
SUPPORT_WIRE_D = 0.0008     # inner support ring wire
SPOKE_N = 3
SPOKE_D = 0.0008
SHAFT_D = 0.0035
HEAD_Z = 0.010              # head mid-plane above the mug bottom

# rotating (MRF) region around the head
MRF_R = 0.013
MRF_Z0 = 0.006
MRF_Z1 = 0.014

# --- mug (inner fluid boundary) ---------------------------------------------
R2_IDEAL = 0.035            # README value, used by the Level 1 idealized annulus
MUG_R_BOTTOM = 0.030        # ESTIMATED
MUG_R_RIM = 0.035
MUG_H = 0.105
LIQUID_H = 0.040


def mug_radius(z):
    """Linear taper from the bottom to the rim."""
    return MUG_R_BOTTOM + (MUG_R_RIM - MUG_R_BOTTOM) * z / MUG_H


def couette(r, omega=OMEGA, r1=R1, r2=R2_IDEAL):
    """Exact laminar azimuthal velocity between a rotating inner and a fixed outer cylinder."""
    a = -omega * r1**2 / (r2**2 - r1**2)
    b = omega * r1**2 * r2**2 / (r2**2 - r1**2)
    return a * r + b / r


def couette_torque_per_length(omega=OMEGA, mu=MU, r1=R1, r2=R2_IDEAL):
    """Exact torque per unit length on the inner cylinder (magnitude)."""
    return 4 * math.pi * mu * omega * r1**2 * r2**2 / (r2**2 - r1**2)


def reynolds(omega=OMEGA, nu=NU, r1=R1, r2=R2_IDEAL):
    return omega * r1 * (r2 - r1) / nu
