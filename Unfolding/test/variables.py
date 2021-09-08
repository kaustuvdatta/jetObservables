#!/usr/bin/env python
import numpy as np
from collections import OrderedDict
import array
from array import array
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

nSubVariables[ 'Jet_tau21' ] = {
            #'bins' : [20],
            'bins' : np.arange( 0.0, 1.1, 0.1 ),
            #'bins' : np.concatenate( ( np.array([ 0, 0.2 ]), np.arange( 0.3, 0.9, 0.06 ), np.array([ 1.]) ) ),
            'label' : 'Leading AK8 jet #tau_{21}',
            'alignLeg' : 'right'
            }

nSubVariables_WSel = OrderedDict()

nSubVariables_WSel[ '_tau_0p5_1' ] = {
            'bins' : array('d',[0., 0.2, 0.28, 0.36, 0.44, 0.52, 0.6, 0.68, 0.76, 0.9]),
            'label' : 'Leading AK8 jet #tau_{1}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_1_1' ] = {
            'bins' : array('d',[0., 0.1, 0.2, 0.26, 0.32, 0.38, 0.44, 0.5, 0.6]),
            'label' : 'Leading AK8 jet #tau_{1}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_2_1' ] = {
            'bins' : array('d',[0., 0.04, 0.08, 0.12, 0.16, 0.2, 0.26, 0.32, 0.4, 0.56]),
            'label' : 'Leading AK8 jet #tau_{1}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_0p5_2' ] = {
            'bins' : array('d',[0., 0.14, 0.2, 0.26, 0.32, 0.4, 0.46,  0.52,  0.58, 0.66]),
            'label' : 'Leading AK8 jet #tau_{2}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_1_2' ] = {
            'bins' : array('d',[0., 0.05, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30, 0.4]),
            'label' : 'Leading AK8 jet #tau_{2}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_2_2' ] = {
            'bins' : array('d',[0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]),
            'label' : 'Leading AK8 jet #tau_{2}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_0p5_3' ] = {
            'bins' : array('d',[0., 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.46, 0.52, 0.6]),
            'label' : 'Leading AK8 jet #tau_{3}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_1_3' ] = {
            'bins' : array('d',[0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.27]),
            'label' : 'Leading AK8 jet #tau_{3}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_2_3' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.04, 0.08, 0.12, 0.16]),
            'label' : 'Leading AK8 jet #tau_{3}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_0p5_4' ] = {
            'bins' : array('d',[0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.24]),
            'label' : 'Leading AK8 jet #tau_{4}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_1_4' ] = {
            'bins' :  array('d',[0., 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25]),
            'label' : 'Leading AK8 jet #tau_{4}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau_2_4' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.064]),
            'label' : 'Leading AK8 jet #tau_{4}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau21' ] = {
            'bins' : array('d',[0., 0.1, 0.18, 0.26, 0.34, 0.42, 0.5, 0.58, 0.66, 0.74, 0.87, 1.]),
            'label' : 'Leading AK8 jet #tau_{2,1}^{(1)}  (WTA #k_{T})',
            'alignLeg' : 'right'
            }
nSubVariables_WSel[ '_tau32' ] = {
            'bins' : array('d',[0., 0.2, 0.3, 0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88, 1.]),
            'label' : 'Leading AK8 jet #tau_{3,2}^{(1)}  (WTA #k_{T})',
            'alignLeg' : 'right'
            }

nSubVariables_topSel = OrderedDict()

nSubVariables_topSel[ '_tau_0p5_1' ] = {
            'bins' : array('d',[0., 0.2, 0.26, 0.32, 0.38, 0.44,  0.5, 0.56, 0.62, 0.68, 0.74, 0.8, 0.9]),
            'label' : 'Leading AK8 jet #tau_{1}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_1_1' ] = {
            'bins' : array('d',[0., 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.46, 0.54, 0.6, 0.7]), 
            'label' : 'Leading AK8 jet #tau_{1}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_2_1' ] = {
            'bins' : array('d',[0., 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.42, 0.48]),
            'label' : 'Leading AK8 jet #tau_{1}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_0p5_2' ] = {
            'bins' : array('d',[0., 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.46, 0.52, 0.58, 0.64, 0.7, 0.8]),
            'label' : 'Leading AK8 jet #tau_{2}^{(0.5)',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_1_2' ] = {
            'bins' : array('d',[0., 0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48]), 
            'label' : 'Leading AK8 jet #tau_{2}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_2_2' ] = {
            'bins' : array('d',[0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2]),
            'label' : 'Leading AK8 jet #tau_{2}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_0p5_3' ] = {
            'bins' : array('d',[0.0, 0.1, 0.14, 0.18, 0.22, 0.26, 0.3, 0.34, 0.38, 0.44, 0.5, 0.6]), 
            'label' : 'Leading AK8 jet #tau_{3}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_1_3' ] = {
            'bins' : array('d',[0., 0.02, 0.04, 0.06, 0.08, 0.11, 0.14, 0.17, 0.2, 0.24, 0.28]),
            'label' : 'Leading AK8 jet #tau_{3}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_2_3' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.05, 0.06, 0.07, 0.09]), 
            'label' : 'Leading AK8 jet #tau_{3}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_0p5_4' ] = {
            'bins' : array('d',[0.0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.42, 0.5]), 
            'label' : 'Leading AK8 jet #tau_{4}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_1_4' ] = {
            'bins' : array('d',[0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24]),
            'label' : 'Leading AK8 jet #tau_{4}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau_2_4' ] = {
            'bins' : array('d',[0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]),
            'label' : 'Leading AK8 jet #tau_{4}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau21' ] = {
            'bins' : array('d',[0., 0.1, 0.18, 0.26, 0.34, 0.42, 0.5, 0.58, 0.66, 0.74, 0.82, 0.9, 1.]),
            'label' : 'Leading AK8 jet #tau_{2,1}^{(1)} (WTA #k_{T})',
            'alignLeg' : 'right'
            }
nSubVariables_topSel[ '_tau32' ] = {
            'bins' : array('d',[0., 0.12, 0.2, 0.28, 0.36, 0.44, 0.52, 0.6, 0.68, 0.76, 0.88, 1.]),
            'label' : 'Leading AK8 jet #tau_{3,2}^{(1)} (WTA #k_{T})',
            'alignLeg' : 'right'
            }


nSubVariables_WSel_test= OrderedDict()

nSubVariables_WSel_test[ '_tau_0p5_1' ] = {
            'bins' : array('d',[0., 0.9]),
            'label' : 'Leading AK8 jet #tau_{1}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_1_1' ] = {
            'bins' : array('d',[0., 0.6]),
            'label' : 'Leading AK8 jet #tau_{1}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_2_1' ] = {
            'bins' : array('d',[0., 0.56]),
            'label' : 'Leading AK8 jet #tau_{1}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_0p5_2' ] = {
            'bins' : array('d',[0., 0.66]),
            'label' : 'Leading AK8 jet #tau_{2}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_1_2' ] = {
            'bins' : array('d',[0., 0.4]),
            'label' : 'Leading AK8 jet #tau_{2}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_2_2' ] = {
            'bins' : array('d',[0., 0.2]),
            'label' : 'Leading AK8 jet #tau_{2}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_0p5_3' ] = {
            'bins' : array('d',[0., 0.6]),
            'label' : 'Leading AK8 jet #tau_{3}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_1_3' ] = {
            'bins' : array('d',[0., 0.28]),
            'label' : 'Leading AK8 jet #tau_{3}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_2_3' ] = {
            'bins' : array('d',[0., 0.16]),
            'label' : 'Leading AK8 jet #tau_{3}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_0p5_4' ] = {
            'bins' : array('d',[0., 0.24]),
            'label' : 'Leading AK8 jet #tau_{4}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_1_4' ] = {
            'bins' :  array('d',[0., 0.26]),
            'label' : 'Leading AK8 jet #tau_{4}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau_2_4' ] = {
            'bins' : array('d',[0., 0.08]),
            'label' : 'Leading AK8 jet #tau_{4}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau21' ] = {
            'bins' : array('d',[0., 1.]),
            'label' : 'Leading AK8 jet #tau_{2,1}^{(1)}  (WTA #k_{T})',
            'alignLeg' : 'right'
            }
nSubVariables_WSel_test[ '_tau32' ] = {
            'bins' : array('d',[0., 1.]),
            'label' : 'Leading AK8 jet #tau_{3,2}^{(1)}  (WTA #k_{T})',
            'alignLeg' : 'right'
            }


nSubVariables_topSel_test = OrderedDict()

nSubVariables_topSel_test[ '_tau_0p5_1' ] = {
            'bins' : array('d',[0., 0.9]),
            'label' : 'Leading AK8 jet #tau_{1}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_1_1' ] = {
            'bins' : array('d',[0., 0.8]), 
            'label' : 'Leading AK8 jet #tau_{1}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_2_1' ] = {
            'bins' : array('d',[0., 0.48]),
            'label' : 'Leading AK8 jet #tau_{1}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_0p5_2' ] = {
            'bins' : array('d',[0., 0.8]),
            'label' : 'Leading AK8 jet #tau_{2}^{(0.5)',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_1_2' ] = {
            'bins' : array('d',[0., 0.48]), 
            'label' : 'Leading AK8 jet #tau_{2}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_2_2' ] = {
            'bins' : array('d',[0., 0.2]),
            'label' : 'Leading AK8 jet #tau_{2}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_0p5_3' ] = {
            'bins' : array('d',[0.0, 0.6]), 
            'label' : 'Leading AK8 jet #tau_{3}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_1_3' ] = {
            'bins' : array('d',[0., 0.28]),
            'label' : 'Leading AK8 jet #tau_{3}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_2_3' ] = {
            'bins' : array('d',[0., 0.1]), 
            'label' : 'Leading AK8 jet #tau_{3}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_0p5_4' ] = {
            'bins' : array('d',[0.0, 0.5]), 
            'label' : 'Leading AK8 jet #tau_{4}^{(0.5)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_1_4' ] = {
            'bins' : array('d',[0., 0.24]),
            'label' : 'Leading AK8 jet #tau_{4}^{(1)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau_2_4' ] = {
            'bins' : array('d',[0., 0.08]),
            'label' : 'Leading AK8 jet #tau_{4}^{(2)}',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau21' ] = {
            'bins' : array('d',[0., 1.]),
            'label' : 'Leading AK8 jet #tau_{2,1}^{(1)} (WTA #k_{T})',
            'alignLeg' : 'right'
            }
nSubVariables_topSel_test[ '_tau32' ] = {
            'bins' : array('d',[0., 1.]),
            'label' : 'Leading AK8 jet #tau_{3,2}^{(1)} (WTA #k_{T})',
            'alignLeg' : 'right'
            }

#nSubVariables[ 'Jet_tau21' ] = {
#            #'bins' : [20],
#            'bins' : np.arange( 0.0, 1.1, 0.1 ),
#            #'bins' : np.concatenate( ( np.array([ 0, 0.2 ]), np.arange( 0.3, 0.9, 0.06 ), np.array([ 1.]) ) ),
#            'label' : 'Leading AK8 jet #tau_{21}',
#            'alignLeg' : 'right'
#            }
