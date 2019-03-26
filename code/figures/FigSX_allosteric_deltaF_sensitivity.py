# -*- coding=utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()

## Define the functions to evaluate the derivatives
def deriv_epai(ep_value, c, ka, ki, n):
    allo_ratio = ((ka * (ki + c))/ (ki * (ka + c)))**n
    return allo_ratio / (allo_ratio + np.exp(ep_value))

def deriv_c(ep, c_value, ka, ki, n):
    numer = -n * ((ka * (ki + c_value)) / (ki * (ka + c_value)))**n * (ka - ki)
    denom_a = (ka + c_value) * (ki + c_value)
    denom_b = ((ka * (ki + c_value) / (ki * (ka + c_value))))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

def deriv_ka(ep, c, ka_value, ki, n):
    numer = -c * n * (ka_value * (ki + c) / (ki * (ka_value + c)))**n
    denom_a = ka_value * (ka_value + c)  
    denom_b = ((ka_value * (ki + c)) / (ki * (ka_value + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)
def deriv_ki(ep, c, ka, ki_value, n):
    numer = c * n * ((ka * (ki_value + c)) / (ki_value * (ka + c)))**n
    denom_a = ki_value * (ki_value + c) 
    denom_b = ((ka * (ki_value + c)) / (ki_value * (ka + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)


# Define ranges to be the same as the interactive zone
n_points = 500
ka_range = np.logspace(-2, 4, n_points)
ki_range = np.logspace(-3, 2, n_points)
epAI_range = np.linspace(-15, 15, n_points)
c_range = np.logspace(-3, 3, n_points)

# Define the constants
c = 50# in µM
ka = 139 # in µM
ki = 0.53 # in µM
ep_ai = 4.5 # in kBT
n = 2

# ##############################
# dF_depAI
# ##############################
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), sharex=True, sharey=True)

# Varying C
d_c_0 = deriv_epai(epAI_range, c * 0.01, ka, ki, n)
d_c_1 = deriv_epai(epAI_range, c, ka, ki, n)
d_c_2 = deriv_epai(epAI_range, c * 100, ka, ki, n)
ax[0].plot(epAI_range, d_c_0, lw=2, color=pboc['blue'], label='0.01')
ax[0].plot(epAI_range, d_c_1, lw=2, color=pboc['red'], label='1')
ax[0].plot(epAI_range, d_c_2, lw=2, color=pboc['purple'], label='100')

# Varying Ka
d_ka_0 = deriv_epai(epAI_range, c, ka * 0.001, ki, n)
d_ka_1 = deriv_epai(epAI_range, c, ka, ki, n)
d_ka_2 = deriv_epai(epAI_range, c, ka * 1000, ki, n)
ax[1].plot(epAI_range, d_ka_0, lw=2, color=pboc['blue'], label='0.001')
ax[1].plot(epAI_range, d_ka_1, lw=2, color=pboc['red'], label='1')
ax[1].plot(epAI_range, d_ka_2, lw=2, color=pboc['purple'], label='1000')

# Varying Ki
d_ki_0 = deriv_epai(epAI_range, c, ka, ki * 0.001, n)
d_ki_1 = deriv_epai(epAI_range, c, ka, ki, n)
d_ki_2 = deriv_epai(epAI_range, c, ka, ki * 1000, n)
ax[2].plot(epAI_range, d_ki_0, lw=2, color=pboc['blue'], label='0.001')
ax[2].plot(epAI_range, d_ki_1, lw=2, color=pboc['red'], label='1')
ax[2].plot(epAI_range, d_ki_2, lw=2, color=pboc['purple'], label='1000')

# Apply special formatting, labels, and legends
leg_titles = ['$c^* / c^{(ref)}$', '$K_A^* / K_A^{(ref)}$', '$K_I^* / K_I^{(ref)}$']
titles = ['varying $c^*$', 'varying $K_A^*$', 'varying $K_I^*$']
panels = ['(A)', '(B)', '(C)']
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_ylim([-0.1, 1.1])
    a.set_xlim([-15, 15])
    _leg = a.legend(loc='lower left', fontsize=7.5, handlelength=0.5, title=leg_titles[i])
    _leg.get_title().set_fontsize(7.5)
    a.set_xlabel(r'∆ε$_{AI}^*$ [$k_BT$]', fontsize=10)
    a.set_title(titles[i], fontsize=8, backgroundcolor='#FFEDC0', y=1.08)

ax[0].set_ylabel(r'$\partial$∆$F / \partial$∆ε$_{AI}^*$', fontsize=10)
fig.text(0, 0.9, '(A)')
fig.text(0.39, 0.9, '(B)')
fig.text(0.68, 0.9, '(C)')
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dEpAI.pdf', bbox_inches='tight')


# ##############################
# dF_dKa
# ##############################
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), sharex=True)

# Varying C
d_c_0 = deriv_ka(ep_ai, c * 0.01, ka_range, ki, n)
d_c_1 = deriv_ka(ep_ai, c, ka_range, ki, n)
d_c_2 = deriv_ka(ep_ai, c * 100, ka_range, ki, n)
ax[0].semilogx(ka_range, d_c_0, lw=2, color=pboc['blue'], label='0.01')
ax[0].semilogx(ka_range, d_c_1, lw=2, color=pboc['red'], label='1')
ax[0].semilogx(ka_range, d_c_2, lw=2, color=pboc['purple'], label='100')

# Varying epAI 
d_ep_0 = deriv_ka(-2, c, ka_range, ki, n)
d_ep_1 = deriv_ka(ep_ai, c, ka_range, ki, n)
d_ep_2 = deriv_ka(8, c, ka_range, ki, n)
ax[1].semilogx(ka_range, d_ep_0, lw=2, color=pboc['blue'], label='-2')
ax[1].semilogx(ka_range, d_ep_1, lw=2, color=pboc['red'], label='4.5')
ax[1].semilogx(ka_range, d_ep_2, lw=2, color=pboc['purple'], label='8')

# Varying Ki
d_ki_0 = deriv_ka(ep_ai, c, ka_range, ki * 0.01, n)
d_ki_1 = deriv_ka(ep_ai,c, ka_range, ki,  n)
d_ki_2 = deriv_ka(ep_ai, c, ka_range, ki * 100, n)
ax[2].semilogx(ka_range, d_ki_0, lw=2, color=pboc['blue'], label='0.01')
ax[2].semilogx(ka_range, d_ki_1, lw=2, color=pboc['red'], label='1')
ax[2].semilogx(ka_range, d_ki_2, lw=2, color=pboc['purple'], label='100')

# Apply special formatting, labels, and legends
leg_titles = ['$c^* / c^{(ref)}$', '∆ε$_{AI}^*$ [$k_BT$]', '$K_I^* / K_I^{(ref)}$']
titles = ['varying $c^*$', 'varying ∆ε$_{AI}^*$', 'varying $K_I^*$']
panels = ['(A)', '(B)', '(C)']
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_ylabel(r'$\partial$∆$F / \partial K_A^*$', fontsize=10)
#     a.set_ylim([-0.1, 1.1])
#     a.set_xlim([-15, 15])
    _leg = a.legend(loc='lower right', fontsize=7.5, handlelength=0.5, title=leg_titles[i])
    _leg.get_title().set_fontsize(7.5)
    a.set_xlabel(r'$K_A^*$ [µM]', fontsize=10)
    a.set_title(titles[i], fontsize=8, backgroundcolor='#FFEDC0', y=1.08)

fig.text(0, 0.9, '(A)')
fig.text(0.39, 0.9, '(B)')
fig.text(0.68, 0.9, '(C)')
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dKa.pdf', bbox_inches='tight')

# ##############################
# dF_dKi
# ##############################
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), sharex=True)

# Varying C
d_c_0 = deriv_ki(ep_ai, c * 0.01, ka, ki_range, n)
d_c_1 = deriv_ki(ep_ai, c, ka, ki_range, n)
d_c_2 = deriv_ki(ep_ai, c * 100, ka, ki_range, n)
ax[0].semilogx(ki_range, d_c_0, lw=2, color=pboc['blue'], label='0.01')
ax[0].semilogx(ki_range, d_c_1, lw=2, color=pboc['red'], label='1')
ax[0].semilogx(ki_range, d_c_2, lw=2, color=pboc['purple'], label='100')

# Varying epAI 
d_ep_0 = deriv_ki(-2, c, ka, ki_range, n)
d_ep_1 = deriv_ki(ep_ai, c, ka, ki_range, n)
d_ep_2 = deriv_ki(8, c, ka, ki_range, n)
ax[1].semilogx(ki_range, d_ep_0, lw=2, color=pboc['blue'], label='-2')
ax[1].semilogx(ki_range, d_ep_1, lw=2, color=pboc['red'], label='4.5')
ax[1].semilogx(ki_range, d_ep_2, lw=2, color=pboc['purple'], label='8')

# Varying Ka
d_ka_0 = deriv_ki(ep_ai, c, ka * 0.01, ki_range, n)
d_ka_1 = deriv_ki(ep_ai,c, ka, ki_range,  n)
d_ka_2 = deriv_ki(ep_ai, c, ka * 100, ki_range, n)
ax[2].semilogx(ki_range, d_ka_0, lw=2, color=pboc['blue'], label='0.01')
ax[2].semilogx(ki_range, d_ka_1, lw=2, color=pboc['red'], label='1')
ax[2].semilogx(ki_range, d_ka_2, lw=2, color=pboc['purple'], label='100')

# Apply special formatting, labels, and legends
leg_titles = ['$c^* / c^{(ref)}$', '∆ε$_{AI}^*$ [$k_BT$]', '$K_A^* / K_A^{(ref)}$']
titles = ['varying $c^*$', 'varying ∆ε$_{AI}^*$', 'varying $K_A^*$']
panels = ['(A)', '(B)', '(C)']
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_yscale('log')
    a.set_ylabel(r'$\partial$∆$F / \partial K_I^*$', fontsize=10)
#     a.set_ylim([-0.1, 1.1])
#     a.set_xlim([-15, 15])
    _leg = a.legend(loc='lower left', fontsize=7.5, handlelength=0.5, title=leg_titles[i])
    _leg.get_title().set_fontsize(7.5)
    a.set_xlabel(r'$K_I^*$ [µM]', fontsize=10)
    a.set_title(titles[i], fontsize=8, backgroundcolor='#FFEDC0', y=1.08)


fig.text(0, 0.9, '(A)')
fig.text(0.39, 0.9, '(B)')
fig.text(0.68, 0.9, '(C)')
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dKi.pdf', bbox_inches='tight')

    
# ##############################
# dF_dc
# ##############################
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), sharex=True)

# Varying epAI 
d_ep_0 = deriv_c(-2, c_range, ka, ki, n)
d_ep_1 = deriv_c(ep_ai, c_range, ka, ki, n)
d_ep_2 = deriv_c(8, c_range, ka, ki, n)
ax[0].semilogx(c_range, d_ep_0, lw=2, color=pboc['blue'], label='-2')
ax[0].semilogx(c_range, d_ep_1, lw=2, color=pboc['red'], label='4.5')
ax[0].semilogx(c_range, d_ep_2, lw=2, color=pboc['purple'], label='8')

# Varying Ka
d_ka_0 = deriv_c(ep_ai, c_range, ka * 0.01, ki, n)
d_ka_1 = deriv_c(ep_ai, c_range, ka, ki,  n)
d_ka_2 = deriv_c(ep_ai, c_range, ka * 100, ki, n)
ax[1].semilogx(c_range, d_ka_0, lw=2, color=pboc['blue'], label='0.01')
ax[1].semilogx(c_range, d_ka_1, lw=2, color=pboc['red'], label='1')
ax[1].semilogx(c_range, d_ka_2, lw=2, color=pboc['purple'], label='100')

# Varying Ki
d_ki_0 = deriv_c(ep_ai, c_range, ka, ki * 0.01, n)
d_ki_1 = deriv_c(ep_ai, c_range, ka, ki, n)
d_ki_2 = deriv_c(ep_ai, c_range, ka, ki * 100, n)
ax[2].semilogx(c_range, d_ki_0, lw=2, color=pboc['blue'], label='0.01')
ax[2].semilogx(c_range, d_ki_1, lw=2, color=pboc['red'], label='1')
ax[2].semilogx(c_range, d_ki_2, lw=2, color=pboc['purple'], label='100')

# Apply special formatting, labels, and legends
leg_titles = ['∆ε$_{AI}^*$ [$k_BT$]', '$K_A^* / K_A^{(ref)}$', '$K_I^* / K_I^{(ref)}$']
titles = ['varying ∆ε$_{AI}^*$', 'varying $K_A^*$', 'varying $K_I^*$']
panels = ['(A)', '(B)', '(C)']
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8) 
    a.set_xlabel(r'$c^*$ [µM]', fontsize=10)
    a.set_ylabel(r'$\partial$∆$F / \partial c^*$', fontsize=10)
    if i != 1:
       _leg = a.legend(loc='lower right', fontsize=7.5, handlelength=0.5, title=leg_titles[i])
       _leg.get_title().set_fontsize(7.5)
    a.set_title(titles[i], fontsize=8, backgroundcolor='#FFEDC0', y=1.08)
_leg = ax[1].legend(loc='lower left', fontsize=7.5, handlelength=0.5, title=leg_titles[i])
_leg.get_title().set_fontsize(7.5)

fig.text(0, 0.9, '(A)')
fig.text(0.39, 0.9, '(B)')
fig.text(0.68, 0.9, '(C)')
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dc.pdf', bbox_inches='tight')
 #
    