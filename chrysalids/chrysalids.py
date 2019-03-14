import sys
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import itertools
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
import glob
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')

def mutant_app():
    bins = 75 
    # Define the inducer concentrations. 
    c_range = np.logspace(-2, 4, 500)
    bohr_range = np.linspace(-20, 20, 500)

    # Define the default reference plot. 
    init_states = dict(R=260, ep_r=-13.9, ka=139, ki=0.52, ep_ai=4.5, n_sites=2,
                       effector_conc=c_range, n_ns=4.6E6)
    ref_arch = mut.thermo.SimpleRepression(**init_states)
    ref_bohr = ref_arch.bohr_parameter()
    ref_foldchange = ref_arch.fold_change()

    # Define the data sources 
    ref_fc_source = ColumnDataSource(dict(c=[c_range], fc=[ref_foldchange],
                                      delta_F=[ref_bohr - ref_bohr], bohr=ref_bohr))
    mut_fc_source = ColumnDataSource(dict(c=[c_range], fc=[ref_foldchange],
                                      delta_F=[ref_bohr-ref_bohr], bohr=ref_bohr))
    mut_point_source = ColumnDataSource(dict(c=c_range[::bins], fc=ref_foldchange[::bins],
                                      delta_F=[ref_bohr-ref_bohr][::bins], bohr=ref_bohr[::bins]))

    # #####################
    # INSTANTIATING FIGURES
    # #####################
    p_foldchange = bokeh.plotting.figure(width=400, height=300, 
                                         x_axis_type='log',
                                         x_axis_label='IPTG [µM]',
                                         y_axis_label='fold-change',
                                         y_range=[-0.1, 1.1], title='Induction Profile')

    p_collapse = bokeh.plotting.figure(width=400, height=300, 
                                         x_axis_label='free energy [kT]',
                                         y_axis_label='fold-change',
                                         y_range=[-0.1, 1.2], 
                                         x_range=[-10, 10], title='Phenotypic Data Collapse')

    p_deltaF = bokeh.plotting.figure(width=400, height=300, 
                                         x_axis_type='log',
                                         x_axis_label='IPTG [µM]',
                                         y_axis_label='∆F [kT]',
                                         y_range=[-10, 10], title='Change in free energy')
   
    # ##############
    # FOLD-CHANGE PLOTS
    # ##############
    p_foldchange.line(x='c', y='fc', source=mut_fc_source, color=pboc['red']) 
    p_foldchange.circle(x='c', y='fc', source=mut_point_source, color=pboc['red'],
    size=10) 
    p_foldchange.line(x='c', y='fc', color='black', source=ref_fc_source)

    # ####################
    # BOHR PARAMETER PLOTS 
    # ####################
    p_collapse.circle(x='bohr', y='fc', color=pboc['red'], source=mut_point_source,
    size=10)
    p_collapse.line(x=bohr_range, y=(1 + np.exp(-bohr_range))**-1, color='black')

    # ####################
    # DELTA F PLOTS
    # ###################
    p_deltaF.line(x='c', y='delta_F', source=mut_fc_source, color=pboc['red'])
    p_deltaF.circle(x='c', y='delta_F', source=mut_point_source, color=pboc['red'],
                    size=10)
    p_deltaF.line(x=c_range, y=ref_bohr - ref_bohr, color='black')

    # ####################
    # CONTROL
    # ####################
    # Mutants
    reset = Button(label='Reset to Reference', button_type='success')
    repressors = Slider(title='Repressors per cell', 
                        value=260, start=1, end=2000, step=1)
    binding_energy = Slider(title='DNA binding energy [kT]',
                            value=-13.9, start=-20, end=-3, step=0.1)
    Ka = Slider(title='Ka [µM]',
                            value=139, start=1, end=1000, step=1)
    Ki = Slider(title='Ki [µM]',
                            value=0.53, start=0.1, end=50, step=0.1)
    epAI = Slider(title='allosteric energy differece [kT]',
                            value=4.5, start=-10, end=10, step=0.1)
    Nsites = Slider(title='Number of allosteric sites', 
                     value=2, start=0, end=10, step=1)
    # Reference
    ref_reset = Button(label='Reset to Default', button_type='success')
    ref_repressors = Slider(title='Repressors per cell', 
                        value=260, start=1, end=2000, step=1)
    ref_binding_energy = Slider(title='DNA binding energy [kT]',
                            value=-13.9, start=-20, end=-3, step=0.1)
    ref_Ka = Slider(title='Ka [µM]',
                            value=139, start=1, end=1000, step=1)
    ref_Ki = Slider(title='Ki [µM]',
                            value=0.53, start=0.1, end=50, step=0.1)
    ref_epAI = Slider(title='allosteric energy differece [kT]',
                            value=4.5, start=-10, end=10, step=0.1)
    ref_Nsites = Slider(title='Number of allosteric sites', 
                     value=2, start=0, end=10, step=1)



    # Define the callbacks. 
    def _mut_reset():
        repressors.value = ref_repressors.value 
        binding_energy.value = ref_binding_energy.value
        Ka.value = ref_Ka.value
        Ki.value = ref_Ki.value
        epAI.value = ref_epAI.value
        Nsites = ref_Nsites.value

    def _ref_reset():
        ref_repressors.value = init_states['R']
        ref_binding_energy.value = init_states['ep_r']
        ref_Ka.value = init_states['ka']
        ref_Ki.value = init_states['ki']
        ref_epAI.value = init_states['ep_ai']
        ref_Nsites = init_states['n_sites']
  
    def _instantiate_mut_arch():
        return mut.thermo.SimpleRepression(R=repressors.value, 
                                    ep_r=binding_energy.value, 
                                    ka=Ka.value,  ki=Ki.value, ep_ai=epAI.value, 
                                    n_sites=Nsites.value, n_ns=4.6E6,
                                    effector_conc=c_range)
    def _instantiate_ref_arch():
        return mut.thermo.SimpleRepression(R=ref_repressors.value, 
                                    ep_r=ref_binding_energy.value, 
                                    ka=ref_Ka.value, ki=ref_Ki.value, 
                                    ep_ai=ref_epAI.value, 
                                    n_sites=ref_Nsites.value, n_ns=4.6E6,
                                    effector_conc=c_range)
   
    # Define the master callback.
    def update_mut_canvas(ref_bohr):
        arch = _instantiate_mut_arch()
        fc = arch.fold_change()
        bohr = arch.bohr_parameter()
        delF = ref_bohr - bohr
        mut_fc_source.data = dict(c=c_range, fc=fc, delta_F=delF, bohr=bohr)
        mut_point_source.data = dict(c=c_range[::bins], fc=fc[::bins], 
                                    delta_F=delF[::bins], bohr=bohr[::bins])
        
    def update_ref_canvas():
        arch = _instantiate_ref_arch()
        fc = arch.fold_change()
        ref_bohr = arch.bohr_parameter()
        delF = ref_bohr - ref_bohr
        ref_fc_source.data = dict(c=c_range, fc=fc, delta_F=delF, bohr=ref_bohr)
        return ref_bohr
    
    def update_canvas():
        ref_bohr = update_ref_canvas()
        update_mut_canvas(ref_bohr)

    # Instantiate the sliders
    mut_controls = [reset, repressors, binding_energy, Ka, Ki, epAI, Nsites]
    ref_controls = [ref_reset, ref_repressors, ref_binding_energy, ref_Ka, 
                    ref_Ki, ref_epAI, ref_Nsites]
    for c in mut_controls[1:]:
        c.on_change('value', lambda attr, old, new: update_canvas())
    reset.on_click(_mut_reset)
    ref_reset.on_click(_ref_reset)

    for c in ref_controls[1:]:
        c.on_change('value', lambda attr, old, new: update_canvas())
   
    # Assemble the layout options
    mut_inputs = widgetbox(mut_controls, sizing_mode='scale_width')
    ref_inputs = widgetbox(ref_controls, sizing_mode='scale_width')
    lay = layout([[mut_inputs, ref_inputs], [p_foldchange, p_collapse, p_deltaF]])
    tab = bokeh.models.widgets.Panel(child=lay, title='Theoretical Model')
    return tab        




