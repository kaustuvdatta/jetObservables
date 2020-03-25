import FWCore.ParameterSet.Config as cms
process = cms.Process('NANO')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),
)
process.source.fileNames = [

    #  2017
    # data
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/STtWantitop5finclusiveDecays13TeV-powheg-pythia8TuneCUETP8M2T4/190321_164541/0000/nano102x_on_mini94x_2016_mc_NANO_8.root'
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/WJetsToQQHT400to600qc193jTuneCP513TeV-madgraphMLM-pythia8RunIIFall17MiniAODv2-PU2017/190414_201204/0000/nano102x_on_mini94x_2017_mc_NANO_170.root'
    'root://xrootd-cms.infn.it//store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190321_164456/0000/nano102x_on_mini94x_2016_mc_NANO_12.root'
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190321_164443/0000/nano102x_on_mini94x_2016_mc_NANO_105.root'
    #'/store/user/algomez/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/QCDPt800to1000TuneCUETP8M113TeVpythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2/190328_223023/0000/nano102x_on_mini94x_2016_mc_NANO_9.root'
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToLNuTuneCUETP8M113TeV-amcatnloFXFX-pythia8RunIISummer16MiniAODv3-PUMoriond1794X/190321_164609/0000/nano102x_on_mini94x_2016_mc_NANO_99.root'
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WJetsToLNuTuneCUETP8M113TeV-amcatnloFXFX-pythia8RunIISummer16MiniAODv3-PUMoriond1794X_ext2-v1/190321_164623/0001/nano102x_on_mini94x_2016_mc_NANO_1690.root'
    #'/store/user/kadatta/SingleMuon/SingleMuon_Run2016H-17Jul2018-v1/190328_151923/0000/nano102x_on_mini94x_2016_data_NANO_97.root'
    #'/store/user/rappocc/JetHT/JetHT_Run2017B-17Nov2017-v1/180525_135243/0000/test_data_94X_NANO_90.root'
    #'/store/user/rappocc/JetHT/JetHT_Run2017D-17Nov2017-v1/180330_195538/0000/test_data_94X_NANO_85.root'
    # mc:
    #'/store/user/rappocc/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180330_193826/0000/test94X_NANO_69.root'
    #'root://cms-xrd-global.cern.ch//store/user/asparker/NanoAODJMARTools-skims/nanoskim-JetsandLepton-94XMC-TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8-1of3-trees.root'
    #'root://cmseos.fnal.gov//store/user/srappocc/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/180330_194407/0000/test94X_NANO_397.root'
    #'file:nano102x_on_mini94x_2016_mc_NANO_98.root',
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/SingleMuon/SingleMuon_Run2016G-17Jul2018-v1/190402_145114/0004/nano102x_on_mini94x_2016_data_NANO_4269.root'
    #'/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/STt-channelantitop4finclusiveDecays13TeV-powhegV2-madspin-pythia8TuneCUETP8M1/190321_164510/0000/nano102x_on_mini94x_2016_mc_NANO_99.root',
    #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/BulkGravToWW_narrow_M-3000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/A6F736EF-1515-E711-8C8A-48FD8E2824DB.root'
    # 2016
    # data
    #'/store/user/srappocc/JetHT/JetHT_Run2016H-07Aug17-v1/180329_210546/0000/test_data_80X_NANO_101.root'
]
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#process.options = cms.untracked.PSet()

process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('jetObservables_nanoskim.root'),
        #fakeNameForCrab =cms.untracked.bool(True)
        )
process.out = cms.EndPath(process.output)
