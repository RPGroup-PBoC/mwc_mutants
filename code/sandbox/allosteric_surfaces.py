# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()

# Define a function for computing the delta bohr. 
def dbohr(ref_params, per_params):
    # Compute the pact for reference and perturbed
    ref_pact = mut.thermo.MWC(effector_conc=ref_params['effector_conc'],
                             ka=ref_params['Ka'], ki=ref_params['Ki'],
                             ep_ai=ref_params['ep_AI'], 
                             n_sites=ref_params['n_sites']).pact()
    per_pact = mut.thermo.MWC(effector_conc=per_params['effector_conc'],
                             ka=per_params['Ka'], ki=per_params['Ki'],
                             ep_ai=per_params['ep_AI'], 
                             n_sites=per_params['n_sites']).pact()
    
    # Compute and return the delta Bohr. 
    delta_bohr = np.log(ref_pact/per_pact)
    return delta_bohr


# Define the reference states. 
Ka_0 = 100 # in µM
Ki_0 = 0.5 # in µM
epAI_0 = 5 # in kBT
c_0 = 50 # in µM
ref_params = {'effector_conc':c_0,
             'Ka':Ka_0, 'Ki':Ki_0,
             'ep_AI': epAI_0,
             'n_sites':2}

# Define the ranges of interest
Ka_range = np.logspace(0, 3, 500)
Ki_range = np.logspace(-3, 0, 500)
epAI_range = np.linspace(-5, 20, 500)
c_range = np.logspace(-2, 4, 500)

# Compute the surfaces. 

## Ka & c
ka, c = np.meshgrid(Ka_range, c_range)
per_params = {'effector_conc':c,
             'Ka':ka, 'Ki':Ki_0,
             'ep_AI':epAI_0,
             'n_sites':2}
ka_c = dbohr(ref_params, per_params)

## Ki & c
ki, c = np.meshgrid(Ki_range, c_range)
per_params = {'effector_conc':c,
             'Ka':Ka_0, 'Ki':ki,
             'ep_AI':epAI_0,
             'n_sites':2}
ki_c = dbohr(ref_params, per_params)

## epAI_c
ep, c = np.meshgrid(epAI_range, c_range)
per_params = {'effector_conc':c,
             'Ka':Ka_0, 'Ki':Ki_0,
              'ep_AI':ep,
             'n_sites':2}
ep_c = dbohr(ref_params, per_params)

## Ka & Ki
ka, ki = np.meshgrid(Ka_range, Ki_range)
per_params = {'effector_conc':c_0,
             'Ka':ka, 'Ki':ki,
             'ep_AI':epAI_0,
             'n_sites':2}
ka_ki = dbohr(ref_params, per_params)

## Ka & ep
ep, ka = np.meshgrid(epAI_range, Ka_range)
per_params = {'effector_conc':c_0,
             'Ka':ka, 'Ki':Ki_0,
              'ep_AI':ep,
             'n_sites':2}
ep_ka = dbohr(ref_params, per_params)

## Ki & ep
ep, ki = np.meshgrid(epAI_range, Ki_range)
per_params = {'effector_conc':c_0,
             'Ka':Ka_0, 'Ki':ki,
              'ep_AI':ep,
             'n_sites':2}
ep_ki = dbohr(ref_params, per_params)



fig, ax = plt.subplots(2, 3, figsize=(6,4))

# Plot the images.
ax[0, 0].imshow(ka_c.T, cmap='bone', vmin=-4, vmax=15)
ax[0, 1].imshow(ki_c.T, cmap='bone', vmin=-4, vmax=15)
ax[0, 2].imshow(ep_c.T, cmap='bone', vmin=-4, vmax=15)
ax[1, 0].imshow(ka_ki.T, cmap='bone', vmin=-4, vmax=15)
ax[1, 1].imshow(ep_ka.T, cmap='bone', vmin=-4, vmax=15)
ax[1, 2].imshow(ep_ki.T, cmap='bone', vmin=-4, vmax=20)

# plot the contours. 
c1 = ax[0, 0].contour(ka_c.T, levels=[-3, 0, 7], colors='w', linewidths=1, linestyles='-')
c2 = ax[0, 1].contour(ki_c.T, levels=[-3, 0, 7], colors='w', linewidths=1, linestyles='-')
c3 = ax[0, 2].contour(ep_c.T, levels=[-3, 0, 7], colors='w', linewidths=1, linestyles='-')
c4 = ax[1, 0].contour(ka_ki.T, levels=[-3, 0, 7], colors='w', linewidths=1, linestyles='-')
c5 = ax[1, 1].contour(ep_ka.T, levels=[-3, 0, 6], colors='w', linewidths=1, linestyles='-')
c6 = ax[1, 2].contour(ep_ki.T, levels=[-3, 0, 7, 15], colors='w', linewidths=1, linestyles='-')

# Label the contours
cs = [c1, c2, c3, c4, c5, c6]
positions = [[(200, 200), (300, 400)],
             [(100, 200), (200,275), (350, 200)],
             [(100, 100), (250, 150), (350, 50)],
             [(200, 400), (300, 200), (450, 50)],
             [(100, 200), (300, 150), (300, 20)],
             [(100, 100), (50, 250), (100, 400), (300, 350)]]
for a, c, p in zip(ax.ravel(), cs, positions):
    cl = a.clabel(c, c.levels, inline=True, fmt='$\,$ %d $k_BT$',
            manual=p, fontsize=8, inline_spacing=3)
    for l in cl:
        l.set_rotation(0)

# Format ticks and ticklabels. 

## c0 / c
for i in range(3):
    ax[0, i].set_xticks(np.arange(0, 500, 150))
    ax[0, i].set_xticklabels(['$10^{2}$', '$10^{1}$',
                      '$10^0$', '$10^{-1}$'])
    
## epAI difference. (lots of hardcoded garbage here.)
_diff = [np.diff([epAI_0, a])[0] for a in epAI_range]
_diff = np.round(_diff, decimals=0)
_diff = np.array(_diff).astype(int)
_diff_target = _diff[::50]
pos = np.arange(0, 500, 50)
for i in range(2):
    ax[1, i+1].set_xticks(pos[::2])
    ax[1, i+1].set_xticklabels(_diff_target[::2])
ax[0, 2].set_yticks(pos[::2])
ax[0, 2].set_yticklabels(_diff_target[::2])

## Transforming Ka/Ka_0
tform = np.round(np.log10(Ka_0 / Ka_range), decimals=2)
logs = [2, 1, 0, -1]
inds = [np.where(tform == ind)[0][0] for ind in logs]
ax[0, 0].set_yticks(inds)
ax[0, 0].set_yticklabels(['$10^{%s}$' %a for a in logs])
ax[1, 0].set_yticks(inds)
ax[1, 0].set_yticklabels(['$10^{%s}$' %a for a in logs])
ax[1, 2].set_yticks(inds)
ax[1, 2].set_yticklabels(['$10^{%s}$' %a for a in logs])

## Transforming Ka/Ka_0
tform = np.round(np.log10(Ki_0 / Ki_range), decimals=2)
logs = [2, 1, 0]
inds = [np.where(tform == ind)[0][0] for ind in logs]
ax[0, 1].set_yticks(inds)
ax[0, 1].set_yticklabels(['$10^{%s}$' %a for a in logs])
ax[1, 0].set_xticks(inds)
ax[1, 0].set_xticklabels(['$10^{%s}$' %a for a in logs])
ax[1, 1].set_yticks(inds)
ax[1, 1].set_yticklabels(['$10^{%s}$' %a for a in logs])

## Varying Ka & c
ax[0, 0].set_ylabel('$K_{A_0} / K_{A}$')
ax[0, 0].set_xlabel('$c_0 / c$')

## Varying Ki & c
ax[0, 1].set_ylabel('$K_{I_0} / K_{I}$')
ax[0, 1].set_xlabel('$c_0 / c$')

## Varying epAI & c
ax[0, 2].set_ylabel(r'$\Delta\varepsilon_{AI_0} - \Delta\varepsilon_{AI}$ [$k_BT$]')
ax[0, 2].set_xlabel('$c_0 / c$')

## Varying Ka & Ki
ax[1, 0].set_ylabel('$K_{A_0} / K_A$')
ax[1, 0].set_xlabel('$K_{I_0} / K_I$')

## Varying Ka & ep_AI
ax[1, 1].set_ylabel('$K_{A_0} / K_A$')
ax[1, 1].set_xlabel(r'$\Delta\varepsilon_{AI_0} - \Delta\varepsilon_{AI}$ [$k_BT$]')

## Varying Ki & ep_AI
ax[1, 2].set_ylabel('$K_{I_0} / K_I$')
ax[1, 2].set_xlabel(r'$\Delta\varepsilon_{AI_0} - \Delta\varepsilon_{AI}$ [$k_BT$]')

# Remove grids.
for a in ax.ravel():
    a.grid(False)
    
for a in ax.ravel():
    a.spines['top'].set_visible(True)
    a.spines['right'].set_visible(True)
       

fig.text(0, 0.95, '(a)', fontsize=8)
fig.text(0.35, 0.95, '(b)', fontsize=8)
fig.text(0.67, 0.95, '(c)', fontsize=8)
fig.text(0, 0.5, '(d)', fontsize=8)
fig.text(0.35, 0.5, '(e)', fontsize=8)
fig.text(0.67, 0.5, '(f)', fontsize=8)
plt.tight_layout()
plt.savefig('allosteric_deltaF.pdf')