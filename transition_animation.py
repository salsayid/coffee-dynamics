import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import os

R1 = 0.0105
R2 = 0.035
Omega1 = 1257
nu = 0.01

N = 200
r = np.linspace(R1, R2, N)

dt = 0.000001
n_steps = 40000
save_every = 400

u = np.zeros(N)
u[0] = Omega1 * R1
u[-1] = 0.0

snapshots = []
times = []

for step in range(n_steps):
    u[0] = Omega1 * R1
    u[-1] = 0.0
    dudr = np.gradient(u, r)
    d2udr2 = np.gradient(dudr, r)
    rhs = nu * (d2udr2 + dudr / r - u / r**2)
    u_new = u.copy()
    u_new[1:-1] = u[1:-1] + dt * rhs[1:-1]
    u_new[0] = Omega1 * R1
    u_new[-1] = 0.0
    u = u_new
    if step % save_every == 0:
        snapshots.append(u.copy())
        times.append(step * dt)

snapshots = np.array(snapshots)
print(f"Solved: {len(snapshots)} frames")

grid_res = 500
x_lin = np.linspace(-R2, R2, grid_res)
y_lin = np.linspace(-R2, R2, grid_res)
X, Y = np.meshgrid(x_lin, y_lin)
R_cart = np.sqrt(X**2 + Y**2)
Theta_cart = np.arctan2(Y, X)

def get_speed_field(u_profile, time_val):
    speed = np.interp(R_cart, r, u_profile, left=np.nan, right=np.nan)


    omega_field = np.interp(R_cart, r, u_profile / r,
                            left=np.nan, right=np.nan)

    phase = Theta_cart - omega_field * time_val * 0.3

    spiral = 0.18 * np.sin(6 * phase)
    speed_visual = speed * (1 + spiral)

    mask = (R_cart < R1) | (R_cart > R2)
    speed_visual[mask] = np.nan
    speed[mask] = np.nan

    return speed_visual, speed

os.makedirs('figures', exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 8), facecolor='black')
ax.set_facecolor('black')
ax.set_aspect('equal')
ax.axis('off')

norm = Normalize(vmin=0, vmax=Omega1 * R1)
cmap = plt.cm.inferno

speed_vis, speed_raw = get_speed_field(snapshots[0], times[0])

im = ax.imshow(speed_vis,
               extent=[-R2, R2, -R2, R2],
               origin='lower',
               cmap='inferno',
               norm=norm,
               interpolation='bilinear',
               zorder=1)

mug_ring = plt.Circle((0, 0), R2, color='#aaaaaa', fill=False,
                       linewidth=2.5, zorder=4)
ax.add_patch(mug_ring)

frother = plt.Circle((0, 0), R1, color='white', fill=True, zorder=5)
ax.add_patch(frother)

outer_mask = patches.Annulus((0, 0), R2 + 0.02, 0.02,
                               color='black', zorder=3)
ax.add_patch(outer_mask)

sm = ScalarMappable(cmap='inferno', norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label('$u_\\theta$ [m/s]', color='white', fontsize=12)
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
cbar.outline.set_edgecolor('white')

time_text = ax.text(0.05, 0.95, 't = 0.000 s',
                    transform=ax.transAxes,
                    color='white', fontsize=13,
                    fontfamily='monospace',
                    verticalalignment='top', zorder=10)

ax.text(0.5, 0.02, 'my mug',
        transform=ax.transAxes,
        color='#888888', fontsize=9,
        ha='center', zorder=10)

ax.text(0, 0, 'frother', color='black', fontsize=7,
        ha='center', va='center', zorder=6, fontweight='bold')

def update(frame_idx):
    snap = snapshots[frame_idx]
    t = times[frame_idx]
    speed_vis, _ = get_speed_field(snap, t)
    im.set_data(speed_vis)
    time_text.set_text(f't = {t:.3f} s')
    return [im, time_text]

ani = animation.FuncAnimation(fig, update,
                               frames=len(snapshots),
                               interval=60,
                               blit=False)

writer = animation.PillowWriter(fps=15)
ani.save('figures/transition_animation.gif', writer=writer, dpi=150)
plt.close()
print("Saved: figures/transition_animation.gif")