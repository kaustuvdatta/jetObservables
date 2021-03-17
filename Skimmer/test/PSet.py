import FWCore.ParameterSet.Config as cms
process = cms.Process('NANO')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),
)
process.source.fileNames = [
    'root://cms-xrd-global.cern.ch///store/user/kadatta/PFNano/106X_v1/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v6-v2/210120_093555/0000/nano106X_on_mini106X_2017_mc_NANO_1.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/QCDPt-15to7000TuneCUETHS1Flat13TeVherwigppRunIISummer16MiniAODv3-PUMoriond1794X/200504_080532/0000/nano102x_on_mini94x_2016_mc_NANO_67.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/QCDPt-15to7000TuneCUETP8M1FlatP613TeVpythia8RunIISummer16MiniAODv3-PUMoriond1794X/190717_164805/0000/nano102x_on_mini94x_2016_mc_NANO_91.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190321_164456/0000/nano102x_on_mini94x_2016_mc_NANO_12.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190321_164443/0000/nano102x_on_mini94x_2016_mc_NANO_534.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190321_164443/0000/nano102x_on_mini94x_2016_mc_NANO_7.root',
    #'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/SingleMuon/SingleMuon_Run2016E-17Jul2018-v1/190429_210147/0000/nano102x_on_mini94x_2016_data_NANO_393.root',
    #'root://xrootd-cms.infn.it//store/user/algomez/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/QCDPt300to470TuneCUETP8M113TeVpythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190328_222907/0000/nano102x_on_mini94x_2016_mc_NANO_12.root',
]
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#process.options = cms.untracked.PSet()

process.output = cms.OutputModule("PoolOutputModule",
        #fileName = cms.untracked.string('jetObservables_histograms.root'),
        fileName = cms.untracked.string('jetObservables_nanoskim.root'),
        fakeNameForCrab =cms.untracked.bool(True)
        )
process.out = cms.EndPath(process.output)
