# -*- coding:utf-8 -*-
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('mut')
mut.viz.plotting_style()

# Load the data and restrict to the double mutants. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class']=='DBL']

# Load the individual mutant summaries. 
kaki_samples = pd.read_csv('../../data/csv/KaKi_only_samples.csv')
kaki_epAI_samples = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
epRA_samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')

# Load the predictions for the double mutants
DBL_pred = pd.read_csv('../../data/csv/DBL_mutant_predictions.csv')

# Define the mutants of interest
DNA_muts = ['Y20I', 'Q21M', 'Q21A']
IND_muts = ['Q294K', 'F164T', 'Q294V']

# Instantiate the figure.
fig, ax = plt.subplots(8, 3, figsize=(3.42, 8.5))
n_draws = int(1E5)
# ###################################
# DATA
# ###################################
for i, dna in enumerate(DNA_muts):
    # Load the DNA binding energy for the given DNA mutant
    epRA_samps = epRA_samples[(epRA_samples['mutant']==dna) & (epRA_samples['operator']=='O2') &\
                             (epRA_samples['repressors']==260)] 
    epRA_draws = epRA_samps['ep_RA'].sample(n_draws, replace=True)                       
    epRA_mode = epRA_samps.iloc[np.argmax(epRA_samps['lp__'].values)]['ep_RA']
    for j, ind in enumerate(IND_muts):
        
        # Define the mutant and isolate the data. 
        mutant = f'{dna}-{ind}' 
        mut_data = data[data['mutant']==mutant]
        
        # Define the IPTG concentration range. 
        c_range = mut_data['IPTGuM']
        
        # Determine which samples should be used based on inducer mutant
        if ind == 'Q294K':
            _kaki_samples = kaki_epAI_samples[(kaki_epAI_samples['mutant']==ind) & 
                                              (kaki_epAI_samples['operator']=='O2')]
            _kaki_draws = _kaki_samples.sample(n_draws, replace=True)
        else:
            _kaki_samples = kaki_samples[(kaki_samples['mutant']==ind) & (kaki_samples['operator']=='O2')]
            _kaki_samples['ep_AI'] = constants['ep_AI']
            _kaki_draws = _kaki_samples.sample(n_draws, replace=True)
        
        # Isolate modes and samples for each parameter. 
        _kaki_mode = _kaki_samples.iloc[np.argmax(_kaki_samples['lp__'].values)]
        Ka_mode = _kaki_mode['Ka']
        Ki_mode = _kaki_mode['Ki']
        epAI_mode = _kaki_mode['ep_AI'] 
         
        # Compute the predicted Bohr parameter. 
        bohr_pred = mut.thermo.SimpleRepression(R=260, ep_r=epRA_mode, 
                                                ka=Ka_mode, ki=Ki_mode, ep_ai=epAI_mode,
                                               n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                               effector_conc=c_range).bohr_parameter() 
        
        # Compute the measured Bohr Parameter.  
        bohr_meas = -np.log(-1  + (1 / mut_data['fold_change']))            

        # Compute the wt bohr parameter
        bohr_wt = mut.thermo.SimpleRepression(R=260, ep_r=constants['O2'],
                                ka=constants['Ka'], ki=constants['Ki'], 
                                ep_ai=constants['ep_AI'], n_sites=constants['n_sites'], 
                                n_ns=constants['Nns'], effector_conc=c_range).bohr_parameter()
        
        # Compute the delta Bohr paramter.  
        delta_bohr_meas = bohr_meas - bohr_wt
        delta_bohr_pred = bohr_pred - bohr_wt
        
        # Compute the HPD for the predicted bohr. 
        bohr_cred_region = np.zeros((2, len(c_range)))
        dbohr_cred_region = np.zeros((2, len(c_range)))
        for k, c in enumerate(c_range):
            _bohr_pred = mut.thermo.SimpleRepression(R=260, ep_r=epRA_draws.values,
                                            ka=_kaki_draws['Ka'].values, ki=_kaki_draws['Ki'].values, 
                                            ep_ai=_kaki_draws['ep_AI'].values, n_sites=constants['n_sites'], 
                                            n_ns=constants['Nns'], effector_conc=c).bohr_parameter()
            bohr_cred_region[:, k] = mut.stats.compute_hpd(_bohr_pred, 0.95)
            
            # Compute the delta bohr credible region. Values could be negative so this is tricky.   
            for z in range(2):
                dbohr_cred_region[z, k] = np.diff([bohr_wt.iloc[k], bohr_cred_region[z, k]])
                     
             
        # Generate a dataframe for easier plotting. 
        _df = pd.DataFrame([])
        _df['dbohr_meas'] = delta_bohr_meas
        _df['dbohr_pred'] = delta_bohr_pred
        _df['dbohr_pred_min'] = dbohr_cred_region[0, :]
        _df['dbohr_pred_max'] = dbohr_cred_region[1, :]
        _df['IPTGuM'] = c_range
        _df['bohr_pred'] = bohr_pred 
        _df['bohr_pred_min'] = bohr_cred_region[0, :]
        _df['bohr_pred_max'] = bohr_cred_region[1, :]
        _df['fold_change'] = mut_data['fold_change']
        
        # Generate summary points. 
        for g, d in _df.groupby('IPTGuM'):
            # Plot the collapse curve and data. 
            ax[i, j].errorbar(d['bohr_pred'].mean(), d['fold_change'].mean(), 
                              np.std(d['fold_change'])**2 / np.sqrt(len(d)),
                              fmt='.', lw=0.75, ms=2, color=colors[mutant]) 
            ax[i, j].hlines(d['fold_change'].mean(), d['bohr_pred_min'], d['bohr_pred_max'],
                            color=colors[mutant], lw=0.75)
            
            # Plot the predicted and measured delta borh. 
            ax[i+5, j].errorbar(d['dbohr_pred'].mean(), d['dbohr_meas'].mean(), 
                                np.std(d['dbohr_meas'])**2/np.sqrt(len(d)), 
                                lw=0.5, fmt='.', color=colors[mutant], ms=2)
            ax[i+5, j].hlines(d['dbohr_meas'].mean(), d['dbohr_pred_min'], d['dbohr_pred_max'],
                            color=colors[mutant], lw=0.75)    
            
# ###################################
# MASTER CURVES
# ###################################
bohr_range = np.linspace(-10, 10, 200)
dbohr_range = np.linspace(-10, 10, 200)
for i in range(3):
    for j in range(3):
         # Bohr collapse curve       
        ax[i, j].plot(bohr_range, 1 / (1 + np.exp(-bohr_range)), 'k-', lw=0.75)
                      
        # Delta Bohr plot# Instantiate the figure.
        ax[i+5, j].plot(dbohr_range, dbohr_range, 'k-', lw=0.75)

# ####################################
# FORMATTING 
# ###################################

# Format the labeling
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6.5)
    a.yaxis.set_tick_params(labelsize=6.5)

for i in range(3):    
    # Clear a row in the middle for modification in illustrator
    ax[3, i].axis('off')
    ax[4, i].axis('off')
    
    # Add axis labels. 
    ax[i, 0].set_ylabel('fold-change', fontsize=6.5)
    ax[2, i].set_xlabel('Bohr parameter\n[$k_BT$]', fontsize=6.5)
    ax[-1, i].set_xlabel('$\Delta F_\mathrm{pred.}$ [$k_BT$]', fontsize=6.5) 
    ax[i+5, 0].set_ylabel('$\Delta F_\mathrm{meas.}$ [$k_BT$]', fontsize=6.5)
    
    # Set limits and turn off unnecessary labels. 
    for j in range(3):
        ax[i, j].set_xlim([-10, 10])
        ax[i, j].set_xticks([-8, 0, 8])
        ax[i, j].set_ylim([-0.05, 1.2])
        ax[i+5, j].set_xlim([-10, 10]) 
        ax[i+5, j].set_xticks([-8, 0, 8])
        ax[i+5, j].set_yticks([-8, 0, 8])
        ax[i+5, j].set_ylim([-10, 10])
        if (i<2):
            ax[i, j].xaxis.set_ticklabels([])
            ax[i+5, j].xaxis.set_ticklabels([])
        if (j > 0):
            ax[i, j].yaxis.set_ticklabels([])
            ax[i+5,j].yaxis.set_ticklabels([])
            
    # Add column and row labels.
    ax[0, i].set_title(IND_muts[i], fontsize=6.5, y=1.05, backgroundcolor='#FFEDC0')
    ax[i, 0].text(-0.8, 0.5, DNA_muts[i], fontsize=6.5, rotation='vertical', backgroundcolor='#FFEDC0', transform=ax[i, 0].transAxes)
    ax[5,i].set_title(IND_muts[i], fontsize=6.5, y=1.05, backgroundcolor='#FFEDC0')
    ax[i+5, 0].text(-0.8, 0.5, DNA_muts[i], fontsize=6.5, rotation='vertical', backgroundcolor='#FFEDC0', transform=ax[i+5, 0].transAxes)   


# ########################################
# SAVING AND POSITIONING
# ########################################
plt.subplots_adjust(hspace=0.08, wspace=0.08)            
plt.savefig('Fig5_doubles.svg')             
