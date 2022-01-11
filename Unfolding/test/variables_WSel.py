#!/usr/bin/env python
import numpy as np
from collections import OrderedDict
import array
from array import array
#### BE CAREFUL WITH BINS! for unfolding the bins array will by divided by 2 to have more reco than gen bins.


nSubVariables_WSel = OrderedDict()

nSubVariables_WSel[ 'Jet_tau_0p5_1' ] = {
            'bins' : array('d',[0., 0.2, 0.28, 0.36, 0.44, 0.52, 0.6, 0.68, 0.76, 0.9]),
            'label' : '#tau_{1}^{(0.5)}',
            'alignLeg' : 'left'
            }
nSubVariables_WSel[ 'Jet_tau_1_1' ] = {
            'bins' : array('d',[0., 0.1, 0.2, 0.26, 0.32, 0.38, 0.44, 0.5, 0.6]),
            'label' : '#tau_{1}^{(1)}',
            'alignLeg' : 'left'
            }
nSubVariables_WSel[ 'Jet_tau_2_1' ] = {
            'bins' : array('d',[0., 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.3, 0.36, 0.44]),
            'label' : '#tau_{1}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_0p5_2' ] = {
            'bins' : array('d',[0., 0.14, 0.2, 0.26, 0.32, 0.4, 0.46,  0.54,  0.6, 0.66]),
            'label' : '#tau_{2}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_1_2' ] = {
            'bins' : array('d',[0., 0.05, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30, 0.4]),
            'label' : '#tau_{2}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_2_2' ] = {
            'bins' : array('d',[0., 0.016, 0.032, 0.048, 0.064, 0.1, 0.14, 0.2]),
            'label' : '#tau_{2}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_0p5_3' ] = {
            'bins' : array('d',[0., 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.6]),
            'label' : '#tau_{3}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_1_3' ] = {
            'bins' : array('d',[0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.27]),
            'label' : '#tau_{3}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_2_3' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.15]),
            'label' : '#tau_{3}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_0p5_4' ] = {
            'bins' : array('d',[0.0, 0.08, 0.14, 0.2, 0.26, 0.32, 0.38, 0.44, 0.6]),
            'label' : '#tau_{4}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_1_4' ] = {
            'bins' :  array('d',[0., 0.03, 0.06, 0.09, 0.12, 0.16, 0.24]),
            'label' : '#tau_{4}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau_2_4' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.042, 0.06]),
            'label' : '#tau_{4}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau21' ] = {
            'bins' : array('d',[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.]),
            'label' : '#tau_{2,1}^{(1)}  (OP k_{T})',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ 'Jet_tau32' ] = {
            'bins' : array('d',[0., 0.3, 0.5, 0.62, 0.74, 0.86, 1.]),
            'label' : '#tau_{3,2}^{(1)}  (OP k_{T})',
            'alignLeg' : 'left'
            }
