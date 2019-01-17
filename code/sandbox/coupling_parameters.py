# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()



pact_ref = 0.5

pact_pred = np.linspace(0, 1, 100) 
pact_true = np.linspace(0.5, 1, 100)


df_pred = np.log(pact_ref / pact_pred) 
df_true = np.log(pact_ref / pact_true)

ddf = df_pred - df_true
plt.plot(df_pred, ddf)










