# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 16:48:03 2026

@author: findl
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def parse_complex(series):
    """parse complex strings from matcont into numpy complex."""
    out = []
    for s in series:
        #regstring go woo
        s = re.sub(r'\s*([+-])\s*', r'\1', str(s).strip())  # remove spaces around +/-
        s = s.replace('i', 'j')
        try:
            out.append(complex(s))
        except Exception:
            out.append(np.nan + 0j)
    return np.array(out)



eq   = pd.read_csv("Matcont_output_1.csv", header=None)
x_df = pd.read_csv("x.csv",  header=None)
f_df = pd.read_csv("f.csv",  header=None)

# row 0  = x 
# row 1  = y 
# row 2  = gamma_x 
# row 14 = eigenvalue 1
# row 15 = eigenvalue 2

eq_gamma = eq.iloc[2,  1:].values.astype(float)
eq_x1    = eq.iloc[0,  1:].values.astype(float)
eq_y1    = eq.iloc[1,  1:].values.astype(float)

eig1 = parse_complex(eq.iloc[14, 1:])
eig2 = parse_complex(eq.iloc[15, 1:])

# Stable when both eigenvalues have negative real part
eq_stable = (np.real(eig1) < 0) & (np.real(eig2) < 0)


# x.csv row 323 = gamma values for the LC branch
# x.csv rows 0:322 step 2 = x (CI) values at each mesh point
# x.csv rows 1:322 step 2 = y (Lac) values at each mesh point
# f.csv row 41 = non-trivial Floquet multiplier (row 42 is always ~1, trivial)
# LC is stable if |Floquet multiplier| <= 1

lc_gamma = x_df.iloc[323, :].values.astype(float)
x1_rows  = x_df.iloc[0:322:2, :].values.astype(float)
y1_rows  = x_df.iloc[1:322:2, :].values.astype(float)

lc_x_max = x1_rows.max(axis=0)
lc_x_min = x1_rows.min(axis=0)
lc_y_max = y1_rows.max(axis=0)
lc_y_min = y1_rows.min(axis=0)

floquet  = np.abs(parse_complex(f_df.iloc[41, :]))
#there was some minor numerical noise (values of 0.99999 <1) so there where false positives, lower tolerance
# from the phase portraits I figured they'd be stable though, so maybe we should leave it as is?
# hm. also icl the pandas ~ mask is great, this module has been great for learning new python
#lc_stable = floquet < 1
lc_stable = floquet  < 0.999

# The floquet stuff ended up being useless, so instead I might just plot it as a whole




fig1, ax = plt.subplots(figsize=(5, 5))

#scatter plots for equilibria so that matplotlib doesnt just draw a straight line
ax.scatter(eq_gamma[eq_stable],  eq_x1[eq_stable],  c='b', marker='o', s=10, label='Stable eq.')
ax.scatter(eq_gamma[~eq_stable], eq_x1[~eq_stable], c='r', marker='s', s=10, label='Unstable eq.')

#thin lines for limit cycles
for row in x1_rows:
    ax.plot(lc_gamma[lc_stable],  row[lc_stable],  'b-', linewidth=0.4, alpha=0.5)
    ax.plot(lc_gamma[~lc_stable], row[~lc_stable], 'r-', linewidth=0.4, alpha=0.5)
# then empty plots for the legend
ax.plot([], [], 'g-', label='Stable LC')
ax.plot([], [], 'm-', label='Unstable LC')

ax.set_xlim(0, 0.3)
ax.set_xlabel('gamma_x')
ax.set_ylabel('X (CI)')
ax.legend()
plt.tight_layout()
plt.savefig("plots/bifurcation_x_lines.png")
plt.close()


fig2, ay = plt.subplots(figsize=(5, 5))

ay.scatter(eq_gamma[eq_stable],  eq_y1[eq_stable],  c='b', marker='o', s=10, label='Stable eq.')
ay.scatter(eq_gamma[~eq_stable], eq_y1[~eq_stable], c='r', marker='s', s=10, label='Unstable eq.')

for row in y1_rows:
    ay.plot(lc_gamma[lc_stable],  row[lc_stable],  'b-', linewidth=0.4, alpha=0.5)
    ay.plot(lc_gamma[~lc_stable], row[~lc_stable], 'r-', linewidth=0.4, alpha=0.5)
ay.plot([], [], 'b-', label='Stable LC')
ay.plot([], [], 'r-', label='Unstable LC')

ay.set_xlim(0, 0.3)
ay.set_xlabel('gamma_x')
ay.set_ylabel('Y (Lac)')
ay.legend()
plt.tight_layout()
plt.savefig("plots/bifurcation_y_lines.png")
plt.close()

print("done")