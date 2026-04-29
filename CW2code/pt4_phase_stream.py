# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 12:25:52 2026

@author: findl
"""
# Streamplot and phase plane tutorial from https://www.fabriziomusacchio.com/blog/2024-03-17-phase_plane_analysis/
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def model(t, state, alpha, sigma, gamma_x, gamma_y, tau_y):
    x, y = state
    f = (1 + x**2 + alpha * sigma * x**4) / ((1 + x**2 + sigma * x**4) * (1 + y**4))
    dx = f - gamma_x * x
    dy = (f - gamma_y * y) / tau_y
    return [dx, dy]

base = dict(alpha=11, sigma=2, gamma_x=0.105, gamma_y=0.036, tau_y=5)
gamma_y_up = dict(alpha=11, sigma=2, gamma_x=0.105, gamma_y=0.050, tau_y=5)
gamma_y_down = dict(alpha=11, sigma=2, gamma_x=0.105, gamma_y=0.020, tau_y=5)
gamma_x_up = dict(alpha=11, sigma=2, gamma_x=0.150, gamma_y=0.036, tau_y=5)
gamma_x_down = dict(alpha=11, sigma=2, gamma_x=0.060, gamma_y=0.036, tau_y=5)
alpha_down = dict(alpha=5, sigma=2, gamma_x=0.105, gamma_y=0.036, tau_y=5)
alpha_up = dict(alpha=20, sigma=2, gamma_x=0.105, gamma_y=0.036, tau_y=5)
tau_up = dict(alpha=11, sigma=2, gamma_x=0.105, gamma_y=0.036, tau_y=10)
tau_down = dict(alpha=11, sigma=2, gamma_x=0.105, gamma_y=0.036, tau_y=2)
sigma_up = dict(alpha=11, sigma=4, gamma_x=0.105, gamma_y=0.036, tau_y=5)
sigma_down = dict(alpha=11, sigma=0.5, gamma_x=0.105, gamma_y=0.036, tau_y=5)
gamma_y_near = dict(alpha=11, sigma=2, gamma_x=0.1, gamma_y=0.036, tau_y=5)


param_sets = {
    "1. Base (paper values)": base,
    "2. gamma_y = 0.050 (increased)": gamma_y_up,
    "3. gamma_y = 0.020 (decreased)": gamma_y_down,
    #"4. gamma_x = 0.150 (increased)": gamma_x_up,
    #"5. gamma_x = 0.060 (decreased)": gamma_x_down,
   # "6. alpha = 5 (weak activation)": alpha_down,
    #"7. alpha = 20 (strong activation)": alpha_up,
    "4. tau_y = 10 (slow Lac)": tau_up,
    "5. tau_y = 2 (fast Lac)": tau_down,
    #"10. sigma = 4 (high OR2 binding)": sigma_up,
    #"11. sigma = 0.5 (low OR2 binding)": sigma_down,
    "6. gamma_x = 0.100 (boundary)": gamma_y_near,
}


#trajectories
IC      = [0.5, 0.5]
t_span  = [0, 500]
t_eval  = np.linspace(*t_span, 10000)

x_vals  = np.linspace(0.01, 5.0, 400)
y_vals  = np.linspace(0.01, 4.0, 400)
X, Y    = np.meshgrid(x_vals, y_vals)

# plotting
for label, p in param_sets.items():

    # vector field on grid
    dX, dY = model(None, [X, Y], p['alpha'], p['sigma'], p['gamma_x'], p['gamma_y'], p['tau_y'])

    sol = solve_ivp(model, t_span, IC, t_eval=t_eval, args=(p['alpha'], p['sigma'], p['gamma_x'], p['gamma_y'], p['tau_y']),method='RK45')

    fig, ax = plt.subplots()

    # streamplot from fabrizio
    ax.streamplot(X, Y, dX, dY, color="grey", density=0.5, linewidth = 0.75)

    # nullclines can be plotted using contours, https://phaseportrait.github.io/reference/phaseportrait/nullclines/nullclines/
    cntr1 = ax.contour(X, Y, dX, levels=[0], colors='orange',  linewidths=2)
    cntr2 = ax.contour(X, Y, dY, levels=[0], colors='deepskyblue', linewidths=2)

    # legend doesnt work with nulclines, https://github.com/matplotlib/matplotlib/issues/11134
    # and the colour and like cmap and linewidth system is completely different between contour streamplot and plot its so confusing for no reason
    # colours finally look tolerable

    h1,_ = cntr1.legend_elements()
    h2,_ = cntr2.legend_elements()
    

    # trajectory (the commas are to get the manual legend thing working)
    traj, = ax.plot(sol.y[0], sol.y[1], '-m', linewidth=2, label='Trajectory')
    start, = ax.plot(sol.y[0][0], sol.y[1][0], 'mo', markersize=6, label='Start')

    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 4.0)
    ax.set_xlabel('x (CI, dimensionless)')
    ax.set_ylabel('y (Lac, dimensionless)')
    ax.set_title(label)
    ax.legend([h1[0], h2[0], traj, start], ['dx/dt = 0 (CI nullcline)', 'dy/dt = 0 (Lac nullcline)', 'Trajectory', 'Start'])
    plt.tight_layout()

    plt.savefig(f"phase_plots/phase_stream_{label}.png")
    plt.close()

print("done")

nrows = 2
ncols = 3
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 10))
axes = axes.flatten()

for ax, (label, p) in zip(axes, param_sets.items()):

    # vector field
    dX, dY = model(None, [X, Y],p['alpha'], p['sigma'], p['gamma_x'], p['gamma_y'], p['tau_y'])
    
    # trajectory
    sol = solve_ivp(model, t_span, IC, t_eval=t_eval,args=(p['alpha'], p['sigma'],p['gamma_x'], p['gamma_y'], p['tau_y']),method='RK45')
    
    # streamplot 
    ax.streamplot(X, Y, dX, dY,color="grey",density=0.4,linewidth=0.6)
    
    # nullclines
    ax.contour(X, Y, dX, levels=[0],colors='orange', linewidths=1.5)
    ax.contour(X, Y, dY, levels=[0],colors='deepskyblue', linewidths=1.5)

    # trajectory
    ax.plot(sol.y[0], sol.y[1], '-m', linewidth=1.5)
    ax.plot(sol.y[0][0], sol.y[1][0], 'mo', markersize=4)

    ax.set_title(label, fontsize=16)
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 4.0)
    ax.tick_params(labelsize=8)

# hide unused axes
for ax in axes[len(param_sets):]:
    ax.set_visible(False)

# shared labels
fig.supxlabel('x (CI, dimensionless)', fontsize=20)
fig.supylabel('y (Lac, dimensionless)', fontsize=20)

# manual legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='orange', lw=2, label='dx/dt = 0'),
    Line2D([0], [0], color='deepskyblue', lw=2, label='dy/dt = 0'),
    Line2D([0], [0], color='m', lw=2, label='Trajectory for {x,y} = {0.5,0.5}'),
    Line2D([0], [0], marker='o', color='m', linestyle='None', label='Start')
]

fig.legend(handles=legend_elements, loc='upper right', fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.8])
plt.savefig("phase_plots/grid_pp.png")
plt.close()

print("done")
