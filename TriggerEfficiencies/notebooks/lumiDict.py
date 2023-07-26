#!/usr/bin/env python
from collections import OrderedDict
#### all reported lumis are from brilcalc and are un /ub 
dictSamples = {

    '2016_preVFP' :  {
            
            'lumiTot' : 19497.905,
            'triggerLumis' : {
                'AK8PFJet80' : 27661.64,
                'AK8PFJet140' : 2689.77,
                'AK8PFJet200' : 278.62,
                'AK8PFJet260' : 56.78,
                'AK8PFJet320' : 19.59,
                'AK8PFJet400' : 6.41,
                'AK8PFJet450' : 1.,
                'AK8PFJet500' : 1.,
                }
            },
        '2016' :  {
            
            'lumiTot' : 16810.813,
            'triggerLumis' : {
                'AK8PFJet80' : 41995.50,
                'AK8PFJet140' : 4319.21,
                'AK8PFJet200' : 652.64,
                'AK8PFJet260' : 75.17,
                'AK8PFJet320' : 25.00,
                'AK8PFJet400' : 8.47,
                'AK8PFJet450' : 1.,
                'AK8PFJet500' : 1.,
                }
            },
        '2017' :  {
            
            'lumiTot' : 41527.74,
            'triggerLumis' : {
                'AK8PFJet80' : 16419.91,
                'AK8PFJet140' : 1560.15,
                'AK8PFJet200' : 219.70,
                'AK8PFJet260' : 88.43,
                'AK8PFJet320' : 33.827,
                'AK8PFJet400' : 5.404,
                'AK8PFJet450' : 4.299,
                'AK8PFJet500' : 1.,
                }
            },
        
        '2018' :  {
            
            'lumiTot' : 59817.407,
            'triggerLumis' : {
                'AK8PFJet60' : 27585.68,
                'AK8PFJet80' : 27585.68,
                'AK8PFJet140' : 1268.74,
                'AK8PFJet200' : 295.24,
                'AK8PFJet260' : 127.99,
                'AK8PFJet320' : 48.17,
                'AK8PFJet400' : 16.08,
                'AK8PFJet450' : 8.09,
                'AK8PFJet500' : 1.,
                }
            },
                    
        }
dictSamples = OrderedDict(dictSamples)