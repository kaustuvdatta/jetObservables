#!/usr/bin/env python
import numpy as np
from collections import OrderedDict

#### BE CAREFUL WITH BINS! for unfolding the bins array will by divided by 2 to have more reco than gen bins.

nSubVariables = OrderedDict()
for ijet in [ ('Jet1', 'Outer' ), ('Jet2', 'Central')]:
    nSubVariables[ ijet[0]+'_tau21' ] = {
                'bins' : np.concatenate( ( np.array([ 0, 0.3 ]), np.arange( 0.36, 0.9, 0.06 ), np.array([ 1.]) ) ),
                'label' : ijet[1]+' AK8 jet #tau_{21}',
                'alignLeg' : 'left'
                }
    nSubVariables[ ijet[0]+'_tau32' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.52 ]), np.arange( 0.58, 1., 0.06 ) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{32}',
		'alignLeg' : 'left'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_1' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.16 ]), np.arange( 0.22, 0.71, 0.06 ), np.array([ 1.]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_2' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.14 ]), np.arange( 0.18, 0.46, 0.04 ), np.array([ .7 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_3' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.1 ]), np.arange( 0.14, 0.40, 0.04 ), np.array([ .6 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_4' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.08 ]), np.arange( 0.12, 0.38, 0.04 ), np.array([ .6 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_1' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.48, 0.06 ), np.array([ .7 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_2' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.16, 0.02 ), np.array([ .5 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_3' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.14, 0.02 ), np.array([ .36 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_4' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.12, 0.02 ), np.array([ .3 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_1' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.24, 0.04 ), np.array([ .5 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_2' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.1, 0.02 ), np.array([ .3 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_3' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.05, 0.01 ), np.array([ .2 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_4' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.04, 0.01 ), np.array([ .15 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(2)}',
		'alignLeg' : 'right'
		}

for ijet in [ ('sdJet1', 'Outer SD' ), ('sdJet2', 'Central SD' )]:
    nSubVariables[ ijet[0]+'_tau21' ] = {
                'bins' : np.concatenate( ( np.array([ 0, 0.2 ]), np.arange( 0.3, 0.78, 0.06 ), np.array([ 1.]) ) ),
                'label' : ijet[1]+' AK8 jet #tau_{21}',
                'alignLeg' : 'left'
                }
    nSubVariables[ ijet[0]+'_tau32' ] = {
		'bins' :  np.concatenate( ( np.array([ 0, 0.4 ]), np.arange( 0.46, .94, 0.06 ), np.array ([1.]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{32}',
		'alignLeg' : 'left'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_1' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.64, 0.08 ), np.array([ 1.]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_2' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.48, 0.06 ), np.array([ .7 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_3' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.42, 0.06 ), np.array([ .6 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_0p5_4' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.42, 0.06 ), np.array([ .6 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(0.5)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_1' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.42, 0.06 ), np.array([ .7 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_2' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.20, 0.04 ), np.array([ .5 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_3' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.20, 0.04 ), np.array([ .36 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_1_4' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.16, 0.02 ), np.array([ .3 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(1)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_1' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., .24, .04), np.array([ .5 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{1}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_2' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.08, 0.02 ), np.array([ .3 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{2}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_3' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.06, 0.02 ), np.array([ .2 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{3}^{(2)}',
		'alignLeg' : 'right'
		}
    nSubVariables[ ijet[0]+'_tau_2_4' ] = {
		'bins' :  np.concatenate( ( np.arange( 0., 0.04, 0.01 ), np.array([ .15 ]) ) ),
		'label' : ijet[1]+' AK8 jet #tau_{4}^{(2)}',
		'alignLeg' : 'right'
		}

#nSubVariables[ 'Jet_tau21' ] = {
#            #'bins' : [20],
#            'bins' : np.arange( 0.0, 1.1, 0.1 ),
#            #'bins' : np.concatenate( ( np.array([ 0, 0.2 ]), np.arange( 0.3, 0.9, 0.06 ), np.array([ 1.]) ) ),
#            'label' : 'Leading AK8 jet #tau_{21}',
#            'alignLeg' : 'right'
#            }
