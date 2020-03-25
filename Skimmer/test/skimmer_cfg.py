import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets

process = cms.Process("Demo")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)  )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(['/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/STt-channelantitop4finclusiveDecays13TeV-powhegV2-madspin-pythia8TuneCUETP8M1/190321_164510/0000/nano102x_on_mini94x_2016_mc_NANO_99.root',
            ])
        )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

####### CMS way to save output root files
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("skimOutputFile_Wprime_M-1000_1200_1600_3500.root"),
        closeFileFast = cms.untracked.bool(True)
        )


####### Selecting genparticles without muons
process.GenParticlesNoNu = cms.EDFilter("CandPtrSelector",
        src = cms.InputTag("packedGenParticles"),
        cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
        )

###### CLusteting genJets, and applying pt cut
process.ak8GenJets = ak4GenJets.clone(
#process.ak8GenJets = ak4PFJets.clone(
        #src = "packedPFCandidates",
        src = "GenParticlesNoNu",
        rParam = 0.8,
        jetAlgorithm = 'AntiKt',
        jetPtMin = 200,
        )


###### Calculating Nsubjetiness for GenJet collection
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
process.newNsub0p5 = Njettiness.clone(
        src = cms.InputTag( 'ak8GenJets' ),
        Njets=cms.vuint32(range(1,8)),         # compute 1-, 2-, 3-, 4- subjettiness
        # variables for measure definition :
        #measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
        beta = cms.double(0.5),              # CMS default is 1
        R0 = cms.double( 0.8 ),              # CMS default is jet cone size
        axesDefinition = cms.uint32( 0 ),    # CMS default is 1-pass KT axes (6), KT_Axes (0)
        )
process.newNsub1 = process.newNsub0p5.clone( beta = 1.0 )
process.newNsub2 = process.newNsub0p5.clone( beta = 2.0 )

###### Running analyzer, which creates a tree with variables
process.addingNsub = cms.EDAnalyzer("Skimmer",
        GenjetNsub0p5 = cms.string( 'newNsub0p5' ) ,
        GenjetNsub1 = cms.string( 'newNsub1' ) ,
        GenjetNsub2 = cms.string( 'newNsub2' ) ,
        srcGenJet = cms.InputTag( 'ak8GenJets'),
        )


process.path = cms.Path(
        process.GenParticlesNoNu
        * process.ak8GenJets
        * process.newNsub0p5
        * process.newNsub1
        * process.newNsub2
        * process.addingNsub
         )
