#%%
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

# %% 
# Determine the EC50
# Define the reference states
Ka_0 = 100 # in µM
Ki_0 = 0.5 # in µM
epAI_0 = 5 # in kBT
ec50 = mut.thermo.SimpleRepression(R=100, ep_r=-14, effector_conc=1, ka=Ka_0,
                                   ki=Ki_0, ep_ai=epAI_0, n_sites=2).ec50()
c_0 = np.array([0, ec50, 1E5])
ref_params = {
             'Ka':Ka_0, 'Ki':Ki_0,
             'ep_AI': epAI_0,
             'n_sites':2}

# Define the imporatant ranges. 
Ka_range = np.logspace(-3, 3, 500)
Ki_range = np.logspace(-3, 3, 500)
epAI_range = np.linspace(-5, 15, 500)

# Compute the surfaces
## Ka & Ki
ka, ki, c = np.meshgrid(Ka_range, Ki_range, c_0)
per_params = {'effector_conc':c,
             'Ka':ka, 'Ki':ki,
             'ep_AI':epAI_0,
             'n_sites':2}
ref_params['effector_conc'] = c
ka_ki = dbohr(ref_params, per_params)

## Ka & ep
ep, ka, c = np.meshgrid(epAI_range, Ka_range, c_0)
per_params = {'effector_conc':c,
             'Ka':ka, 'Ki':Ki_0,
              'ep_AI':ep,
             'n_sites':2}
ref_params['effector_conc'] = c
ep_ka = dbohr(ref_params, per_params)

## Ki & ep
ep, ki, c = np.meshgrid(epAI_range, Ki_range, c_0)
per_params = {'effector_conc':c,
             'Ka':Ka_0, 'Ki':ki,
              'ep_AI':ep,
             'n_sites':2}
ref_params['effector_conc'] = c
ep_ki = dbohr(ref_params, per_params)

# %%

# Generate the figure. 
fig, ax = plt.subplots(3, 3, figsize=(6, 6))

# Show the heatmaps
cs  =  []
for i in range(3):
    ax[0, i].imshow(ka_ki[:, :, i].T, cmap='bone', vmin=-8, vmax=20)
    ax[1, i].imshow(ep_ka[:, :, i], cmap='bone', vmin=-8, vmax=20)
    ax[2, i].imshow(ep_ki[:, :, i], cmap='bone', vmin=-8, vmax=20)

    # Generate and label contours. 
    if i == 0:
        levels = [-1, 0, 1]
    else:
        levels = [-3, 0, 5]
    c1 = ax[0, i].contour(ka_ki[:,:, i].T,levels, colors='w', linewidths=1, linestyles='-')
    c2 = ax[1, i].contour(ep_ka[:,:,i], levels, colors='w', linewidths=1, linestyles='-')
    c3 = ax[2, i].contour(ep_ki[:,:,i], levels, colors='w', linewidths=1, linestyles='-')
    cs.append([c1, c2, c3])
    # Add contour labels

# Manually adjust contour labels. 

## EpAI
for i in range(2):
    cl = ax[i + 1,  0].clabel(cs[0][i+1], cs[0][i+1].levels, inline=True, fmt='$\,$ %d $k_BT$',
             manual=[[100, 350], [350, 250]], fontsize=8, inline_spacing=2)
    for l in cl:
        l.set_rotation(0)

## KaKi
cl1 = ax[0,  1].clabel(cs[1][0], cs[1][0].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[150, 440], [180, 300], [150, 100]])
cl2 = ax[0,  2].clabel(cs[2][0], cs[2][0].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[150, 440], [100, 300], [150, 100]])

cl3 = ax[1,  1].clabel(cs[1][1], cs[1][1].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[150, 440], [200, 350], [150, 100]])
cl4 = ax[1,  2].clabel(cs[2][1], cs[2][1].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[150, 440], [100, 300], [350, 300]])

cl5 = ax[2,  1].clabel(cs[1][2], cs[1][2].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[150, 100], [150, 300], [300, 450]])
cl6 = ax[2,  2].clabel(cs[2][2], cs[2][2].levels, inline=True, fmt='$\,$ %d $k_BT$',
        fontsize=8, inline_spacing=2, manual=[[75, 300], [250, 300], [300, 200]])

for l in cl1 + cl2 + cl3 + cl4 + cl5 + cl6:
    l.set_rotation(0)

# Apply appropriate units to the axes. 
## Transforming Ka/Ka_0
tform = np.round(np.log10(Ka_0 / Ka_range), decimals=1)
logs = [5, 4, 3, 2, 1, 0, -1]
inds = [np.where(tform == ind)[0][0] for ind in logs]
for i in range(3):
    ax[0, i].set_yticks(inds)
    ax[0, i].set_yticklabels(['$10^{%s}$' %a for a in logs])
    ax[1, i].set_yticks(inds)
    ax[1, i].set_yticklabels(['$10^{%s}$' %a for a in logs])
# ## Transforming Ki/Ki_0
tform = np.round(np.log10(Ki_0 / Ki_range), decimals=1)
logs = [2, 1, 0, -1, -2, -3]
inds = [np.where(tform == ind)[0][0] for ind in logs]
for i in range(3):
    ax[0, i].set_xticks(inds)
    ax[0, i].set_xticklabels(['$10^{%s}$' %a for a in logs])
    ax[2, i].set_yticks(inds)
    ax[2, i].set_yticklabels(['$10^{%s}$' %a for a in logs])

## Transforming dep_AI
_diff = epAI_0 - epAI_range
_diff = np.round(_diff, decimals=0)
_diff = np.array(_diff).astype(int)
_diff_target = _diff[::50]
pos = np.arange(0, 500, 50)
for i in range(2):
    for j in range(3):
        ax[i+1, j].set_xticks(pos[::3])
        ax[i+1, j].set_xticklabels(_diff_target[::3])

# Format the axes as needed. 
for a in ax.ravel():
    a.grid(False)
    a.spines['top'].set_visible(True)
    a.spines['right'].set_visible(True)
         
# Add appropriate titles and lables. 
ax[0, 0].set_title('$c = 0$', fontsize=8, backgroundcolor="#f1f2f6", y=1.08)
ax[0, 1].set_title(r'$c = [EC_{50}]$', fontsize=8, 
                   backgroundcolor="#f1f2f6", y=1.08)
ax[0, 2].set_title(r'$c \rightarrow \infty$', fontsize=8, 
                    backgroundcolor="#f1f2f6", y=1.08)

for i in range(3):
    ax[0, i].set_xlabel('$K_{I_0} / K_I$')
    ax[0, i].set_ylabel('$K_{A_0} / K_A$')
    ax[1, i].set_ylabel('$K_{A_0} / K_A$')
    ax[2, i].set_ylabel('$K_{I_0} / K_I$')
    ax[1, i].set_xlabel(r'$\Delta\varepsilon_{AI_0} - \Delta\varepsilon_{AI}$ [$k_BT$]')
    ax[2, i].set_xlabel(r'$\Delta\varepsilon_{AI_0} - \Delta\varepsilon_{AI}$ [$k_BT$]')

# Add panel labels. 
fig.text(0, 0.95,  '(a)', fontsize=9)
fig.text(0, 0.65,  '(b)', fontsize=9)
fig.text(0, 0.35,  '(c)', fontsize=9)
plt.tight_layout()
plt.savefig('allosteric_surfaces.pdf', bbox_inches='tight')