# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.personal_style()


    r0_r = np.logspace(-6, 6)
    plt.semilogx(r0_r, np.log(r0_r))
