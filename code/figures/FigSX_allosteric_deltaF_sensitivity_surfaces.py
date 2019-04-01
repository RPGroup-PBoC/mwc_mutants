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

def nice_contour(axis, mat, levels, fmt, tolerance=100, 
    canvas_size=n_points, contour_kwargs={}, clabel_kwargs={}):
    # Draw contours. 
    contours = axis.contour(mat, levels, **contour_kwargs)
    labels = axis.clabel(contours, fmt=fmt, **clabel_kwargs)
    man_pos = []
    for cont, lab in zip(contours.collections, labels):
        x = []
        _lab = lab.get_position()
        for j in range(2):
            if (_lab[j] < tolerance):
                x.append(tolerance)
            elif (_lab[j] > canvas_size - tolerance):
                x.append(canvas_size - tolerance)
            else:
                x.append(_lab[j])
        cont.remove() 
        lab.remove() 
        man_pos.append(x)
    _c = axis.contour(mat, levels=levels, **contour_kwargs)
    _cls = axis.clabel(_c, fmt=fmt, manual=man_pos, **clabel_kwargs)
    return [_c, _cls]

# Define ranges to be the same as the interactive zone
n_points = 500
ka_range = np.logspace(-2, 4, n_points)
ki_range = np.logspace(-3, 2, n_points)
epAI_range = np.linspace(-8, 15, n_points)
c_range = np.logspace(-3, 3, n_points)

# Define the constants
c = 50# in µM
ka = 139 # in µM
ki = 0.53 # in µM
ep_ai = 4.5 # in kBT
n = 2


# ###############################################
# dF / dEpAI
# ###############################################
# Mesh the parameters
_ep_ka, _ka = np.meshgrid(epAI_range, ka_range)
_ep_ki, _ki = np.meshgrid(epAI_range, ki_range)
_ep_c, _c = np.meshgrid(epAI_range, c_range)

ticks = {'epAI':[[0, 173, 347, 500], [-8, 0, 8, 15]],
         'Ka':[[0, 83, 166, 249, 332, 415, 500], 
               ['$10^{-2}$', '$10^{-1}$', '$10^0$',
                '$10^1$', '$10^2$', '$10^3$', '$10^4$']],
        'Ki':[[0, 100, 200, 300, 400, 500], 
               ['$10^{-3}$', '$10^{-2}$', '$10^{-1}$',
                '$10^0$', '$10^1$', '$10^2$']],
        'c': [[0, 83, 166, 249, 332, 415, 500], 
               ['$10^{-3}$', '$10^{-2}$', '$10^{-1}$',
                '$10^0$', '$10^1$', '$10^2$', '$10^3$']]}
# Compute the surfaces
dFdE_ka = deriv_epai(_ep_ka, c, _ka, ki, n)
dFdE_ki = deriv_epai(_ep_ki, c, ka, _ki, n)
dFdE_c = deriv_epai(_ep_c, _c, ka, ki, n)


# Set up the figure canvas
fig, ax = plt.subplots(1, 3, figsize=(7, 2.5))
    
_levels = [0.001, 0.01, 0.1, 0.5, 0.95, 0.999]
ax[0].imshow(dFdE_c, origin='lower', cmap='bone', vmin=0.001, vmax=1.25)
_ = nice_contour(ax[0], dFdE_c, levels=_levels, fmt='%.3f', tolerance=60, 
                        contour_kwargs=dict(colors=['w']), 
                       clabel_kwargs=dict(fontsize=7))
ax[1].imshow(dFdE_ka, origin='lower', cmap='bone', vmin=-0.001, vmax=1.25)
_ = nice_contour(ax[1], dFdE_ka, levels=_levels, fmt='%.3f', tolerance=60, 
                        contour_kwargs=dict(colors=['w']), 
                       clabel_kwargs=dict(fontsize=7))
ax[2].imshow(dFdE_ki, origin='lower', cmap='bone', vmin=-0.001, vmax=1.25)
_ = nice_contour(ax[2], dFdE_ki, levels=_levels, fmt='%.3f', tolerance=100, 
                        contour_kwargs=dict(colors=['w']), 
                       clabel_kwargs=dict(fontsize=7))

# Format the axes
for a in ax:
    a.grid(False)
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    a.set_xlabel('∆ε$_{AI}^*$ [$k_BT$]')
    
    # This is hardcoded and awful and I hate it
    a.set_xticks(ticks['epAI'][0])
    a.set_xticklabels(ticks['epAI'][1])
    
ax[0].set_ylabel('$c^*$ [µM]', fontsize=9)
ax[0].set_yticks(ticks['c'][0])
ax[0].set_yticklabels(ticks['c'][1])

ax[1].set_yticks(ticks['Ka'][0])
ax[1].set_yticklabels(ticks['Ka'][1])
ax[1].set_ylabel('$K_A^*$ [µM]', fontsize=9)

ax[2].set_yticks(ticks['Ki'][0])
ax[2].set_yticklabels(ticks['Ki'][1])
ax[2].set_ylabel('$K_I^*$ [µM]', fontsize=9)

fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.35, 0.95, '(B)', fontsize=8)
fig.text(0.68, 0.95, '(C)', fontsize=8)
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dEpAI_surfaces.pdf', bbox_inches='tight')


# ###############################################
# dF / dKa ; dF / dKi
# ###############################################
_ka_ep, _ep_ka= np.meshgrid(ka_range, epAI_range)
_ka_c, _c_ka = np.meshgrid(ka_range, c_range)
_ka_ki, _ki_ka = np.meshgrid(ka_range, ki_range)
dFdKa_ep = deriv_ka(_ep_ka, c, _ka_ep, ki, n)
dFdKa_c = deriv_ka(ep_ai, _c_ka, _ka_c, ki, n)
dFdKa_ki = deriv_ka(ep_ai, c, _ka_ki, _ki_ka, n)

_ki_ep, _ep_ki= np.meshgrid(ki_range, epAI_range)
_ki_c, _c_ki = np.meshgrid(ki_range, c_range)
_ki_ka, _ka_ki = np.meshgrid(ki_range, ka_range)
dFdKi_ep = deriv_ki(_ep_ki, c, ka, _ki_ep, n)
dFdKi_c = deriv_ki(ep_ai, _c_ki, ka, _ki_c, n)
dFdKi_ka = deriv_ki(ep_ai, c, _ka_ki, _ki_ka, n)

# Set up the figure canvas
fig, ax = plt.subplots(2, 3, figsize=(7, 4.5))

# Plot epAI dependence
ax[0,0].imshow(np.log10(np.abs(dFdKa_ep)), origin='lower', cmap='bone', vmin=-9, vmax=3)
_, _l = nice_contour(ax[0,0], np.log10(np.abs(dFdKa_ep)), fmt='$-10^{%d}$', levels=[-4, -2, 0, 1],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)
for l in _l:
    l.set_rotation(0)

ax[1,0].imshow(np.log10(dFdKi_ep), origin='lower', cmap='bone', vmin=-9, vmax=4.5)
_, _l = nice_contour(ax[1,0], np.log10(dFdKi_ep), fmt='$10^{%d}$', levels=[-4, -2, 0, 1, 2, 3],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)
for l in _l:
    l.set_rotation(0)

for i in range(2):
    ax[i, 0].set_ylabel('∆ε$_{AI}^*$ [$k_BT$]', fontsize=9)
    ax[i, 0].set_yticks(ticks['epAI'][0])
    ax[i, 0].set_yticklabels(ticks['epAI'][1])
    
# Plot c dependence
ax[0,1].imshow(np.log10(np.abs(dFdKa_c)), origin='lower', cmap='bone', vmin=-5, vmax=-0.1)
_, _l = nice_contour(ax[0,1], np.log10(np.abs(dFdKa_c)), fmt='$-10^{%d}$', levels=[-4,-3, -2, -1],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)
for l in _l:
    l.set_rotation(0)

ax[1,1].imshow(np.log10(dFdKi_c), origin='lower', cmap='bone', vmin=-4, vmax=2)
_, _l = nice_contour(ax[1,1], np.log10(dFdKi_c), fmt='$10^{%d}$', levels=[-2, -1,  0, 1],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)
for l in _l:
    l.set_rotation(0)

for i in range(2):
    ax[i, 1].set_ylabel('$c^*$ [µM]', fontsize=9)
    ax[i, 1].set_yticks(ticks['c'][0])
    ax[i, 1].set_yticklabels(ticks['c'][1])

    
# K dependence    
ax[0,2].imshow(np.log10(np.abs(dFdKa_ki)), origin='lower', cmap='bone', vmin=-8, vmax=2.5)
_, _l = nice_contour(ax[0,2], np.log10(np.abs(dFdKa_ki)), fmt='$-10^{%d}$', levels=[-7,-3, 0, 1],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)
for l in _l:
    l.set_rotation(0)

ax[1,2].imshow(np.log10(dFdKi_ka), origin='lower', cmap='bone', vmin=-3, vmax=3.5)
_, _l = nice_contour(ax[1,2], np.log10(dFdKi_ka), fmt='$10^{%d}$', levels=[-3, 0, 1, 2],
                                   contour_kwargs=dict(colors=['w'], linestyles=['-']), 
                                clabel_kwargs=dict(fontsize=8, inline_spacing=2), tolerance=50)   
for l in _l:
    l.set_rotation(0)
   

ax[0, 2].set_ylabel('$K_I^*$ [µM]', fontsize=9)
ax[0, 2].set_yticks(ticks['Ki'][0])
ax[0, 2].set_yticklabels(ticks['Ki'][1])
ax[1, 2].set_ylabel('$K_A^*$ [µM]', fontsize=9)
ax[1, 2].set_yticks(ticks['Ka'][0])
ax[1, 2].set_yticklabels(ticks['Ka'][1])   

# Format the axes
for a in ax.ravel():
    a.grid(False)
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    
for i in range(3):
    ax[0, i].set_xlabel('$K_A^*$ [µM]', fontsize=9)
    ax[0,i].set_xticks(ticks['Ka'][0])
    ax[0,i].set_xticklabels(ticks['Ka'][1])   
for i in range(3):
    ax[1, i].set_xlabel('$K_I^*$ [µM]', fontsize=9)
    ax[1,i].set_xticks(ticks['Ki'][0])
    ax[1,i].set_xticklabels(ticks['Ki'][1])   
 
fig.text(0, 0.99, '(A)', fontsize=9)
fig.text(0.35, 0.99, '(B)', fontsize=9)
fig.text(0.68, 0.99, '(C)', fontsize=9)
fig.text(0, 0.5, '(D)', fontsize=9)
fig.text(0.35, 0.5, '(E)', fontsize=9)
fig.text(0.68, 0.5, '(F)', fontsize=9)
plt.tight_layout()
plt.savefig('../../figures/FigSX_dF_dKaKi_surfaces.pdf', bbox_inches='tight')




# Plot the profile of Ka alone
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.loglog(c_range, np.abs(dFdKa_c[:, 250]))
# ax.vlines(, 1E-5, 1E-1)












