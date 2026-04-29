# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:22:35 2026

@author: findl
"""



import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# scipy integrate as per
# both have the middle bit so that can be shared
def model(t, state, alpha, sigma, gamma_x, gamma_y, tau_y):
    x, y = state
    f = (1 + x**2 + alpha * sigma * x**4) / ((1 + x**2 + sigma * x**4) * (1 + y**4))
    dx = f - gamma_x * x
    dy = (f - gamma_y * y) / tau_y
    return [dx, dy]


def simulate(params, t_span=(0, 500), t_points=5000, ic=(0.5, 0.5)):
    t_eval = np.linspace(*t_span, t_points)
    sol = solve_ivp(
        model,
        t_span,
        ic,
        args=(params['alpha'], params['sigma'], params['gamma_x'], params['gamma_y'], params['tau_y']),
        t_eval=t_eval,
        method='RK45',

    )
    return sol.t, sol.y[0], sol.y[1]


# base params all from the paper
# sigma=2, alpha=11 from table I, gamma_x=0.105, gamma_y=0.036 from fig 2d
# tau_y=5 stated several times in the paper

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
    "6. gamma_x = 0.100 (near boundary, outside)": gamma_y_near,
}



# plot

for label, params in param_sets.items():
    t, x, y = simulate(params)

    fig, ax = plt.subplots()
    
    ax.plot(t, y, label='y (Lac)')
    ax.plot(t, x, label='x (CI)')
    
    ax.set_xlabel('Dimensionless time')
    ax.set_ylabel('Dimensionless concentration')
    ax.set_title(label)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"plots/{label}.png")
    plt.close()

print("done")
nrows = 2
ncols =3
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 5))
axes = axes.flatten()

for ax, (label, params) in zip(axes, param_sets.items()):
    t, x, y = simulate(params)
    
    ax.plot(t, y, label='y (Lac)')
    ax.plot(t, x, label='x (CI)')
    
    ax.set_title(label, fontsize=12)
    ax.tick_params(labelsize=10)

for ax in axes[len(param_sets):]:
    ax.set_visible(False)

fig.supxlabel('Dimensionless time', fontsize=16)
fig.supylabel('Dimensionless concentration', fontsize=16)

# shared legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=12)

plt.tight_layout(rect=[0, 0, 0.90, 1]) # leave a bit of space for the legend
plt.savefig("plots/grid.png")
plt.close()

print("done")


# looking for peaks adn troughts
def turning_points(t, x, transient_cut=0.5):
    cut_idx = int(len(t) * transient_cut)
    t_ss = t[cut_idx:]
    x_ss = x[cut_idx:]

    dx = np.diff(x_ss)
    signs = np.sign(dx)

    # sign changes
    sign_changes = np.diff(signs)

    peaks = np.where(sign_changes < 0)[0] #python indexing moment
    troughs = np.where(sign_changes > 0)[0] # also .where is magic

    peak_coords = [(t_ss[i], x_ss[i]) for i in peaks]
    trough_coords = [(t_ss[i], x_ss[i]) for i in troughs]

    return peak_coords, trough_coords
for label, params in param_sets.items():
    t, x, y = simulate(params)

    peaks, troughs = turning_points(t, x)

    print(f"\n{label}")
    print("Peaks:", peaks[:5])
    print("Troughs:", troughs[:5])
