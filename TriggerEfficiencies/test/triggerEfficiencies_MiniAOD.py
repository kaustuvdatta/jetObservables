import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)  )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(['/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/50000/C000CF01-A463-7241-BDCE-C3ADA5854799.root',
'/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/50000/8E2BC804-9678-0247-93B3-C80FF4D604B9.root',
'/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/50000/9CDE44B8-48D6-AE4B-BBB9-193C583F765C.root',
'/store/data/Run2017B/JetHT/MINIAOD/UL2017_MiniAODv2-v1/50000/A4D80D52-2D19-2741-BFC7-EFCC61FB76F0.root',
            ])
        )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

####### CMS way to save output root files
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("triggerEfficienciesMiniOAD.root"),
        closeFileFast = cms.untracked.bool(True)
        )

triggerThresholds = [ '60', '80', '140', '200', '260', '320', '400', '450', '500', '550' ]
listOfTriggers = [ 'HLT_AK8PFJet'+str(x) for x in triggerThresholds  ]

###### Running analyzer, which creates a tree with variables
process.TriggerEfficiencies = cms.EDAnalyzer("TriggerEfficiencies",
        cutAK8jetPt= cms.double( 170 ) ,
        baseTrigger = cms.string( 'HLT_AK8PFJet60' ),
        triggerThresholds = cms.vstring( triggerThresholds ),
        listOfTriggers = cms.vstring( listOfTriggers )
        )


process.path = cms.Path(
        process.TriggerEfficiencies
         )
