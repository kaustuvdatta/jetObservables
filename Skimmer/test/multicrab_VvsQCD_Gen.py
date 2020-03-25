##################################################################
########   TO RUN THIS: python multicrab.py -d dataset -v 9_4_X_v2 -s T3_US_FNALLPC
########   In --dataset, if you specify for instance QCD it runs ALL the keys in datasamples that starts with QCD.
########   DO NOT DO: crab submit multicrab.py
##################################################################

from CRABClient.UserUtilities import config
import argparse, sys, os
from httplib import HTTPException
from CRABAPI.RawCommand import crabCommand


#### General crab parameters
config = config()
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'skimmer_cfg.py'
config.Data.inputDBS = 'global'
config.Data.publication = False
#config.Site.whitelist = ['T2_IT_Legnaro', 'T2_US_Caltech', 'T2_TW_NCHC']
#config.Site.blacklist = ['T2_US_Caltech']
#config.Data.ignoreLocality = True
###############################################################

def submit(config):
	try:
		if args.dryrun: crabCommand('submit', '--dryrun', config = config)
		else: crabCommand('submit', config = config)
	except HTTPException, hte:
		print 'Cannot execute commend'
		print hte.headers


#######################################################################################
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--dataset', action='store', default='all', dest='dataset', help='Sample to process. Example: QCD, BTagMu.' )
	parser.add_argument('-v', '--version', action='store', default='9_4_X_v2', dest='version', help='Version. Example: 9_4_X_v1' )
	parser.add_argument('-s', '--storageSite', action='store', default='T3_CH_PSI', dest='storageSite', help='storageSite. Example: T2_CH_CERN or T3_US_FNAL' )
	parser.add_argument('-D', '--dryrun', action='store_true', default=False, help='To run dryrun crab mode.' )
	parser.add_argument('-t', '--testNoSend', action='store_true', default=False, help='To print crab config without send it. Helpful to debug' )

	try: args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	Samples = {}
	###  Samples[Dataset nickname] = [ 'dataset name', splitting number  ]
	Samples['WtoQQ_1'] = ['/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['WtoQQ_2'] = ['/WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['WtoQQ_3'] = ['/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
    
	Samples['QCD_1'] = [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_2'] = [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_3'] = [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_4'] = [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_5'] = [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_6'] = [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_7'] = [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_8'] = ['/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_9'] = ['/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
	Samples['QCD_10'] = ['/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',10000]
    
    
	#Samples['Wqq'] = ['/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM',10000]
					#Samples['ZJetsToQQ_0'] = ['/ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/algomez-ZJetsToQQHT400to600qc194jTuneCP513TeV-madgraphMLM-pythia8RunIIFall17MiniAODv2-PU2017-2632477341b0033d0ee33ee9e5481e57/USER',10000]
	#Samples['ZZto2Q2Nu'] = ['/ZZTo2Q2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', 10000]
	#Samples['WJetsToQQ_0'] = ['/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/algomez-WJetsToQQHT400to600qc193jTuneCP513TeV-madgraphMLM-pythia8RunIIFall17MiniAODv2-PU2017-2632477341b0033d0ee33ee9e5481e57/USER',10000]


#Samples[ 'WprimeToWZToWhadZinv0' ] = [ '/WprimeToWZToWhadZinv_narrow_M-600_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv1' ] = [ '/WprimeToWZToWhadZinv_narrow_M-800_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv2' ] = [ '/WprimeToWZToWhadZinv_narrow_M-1000_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv3' ] = [ '/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv4' ] = [ '/WprimeToWZToWhadZinv_narrow_M-1400_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv5' ] = [ '/WprimeToWZToWhadZinv_narrow_M-1600_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv6' ] = [ '/WprimeToWZToWhadZinv_narrow_M-1800_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv7' ] = [ '/WprimeToWZToWhadZinv_narrow_M-2000_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv8' ] = [ '/WprimeToWZToWhadZinv_narrow_M-2500_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv9' ] = [ '/WprimeToWZToWhadZinv_narrow_M-3000_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'WprimeToWZToWhadZinv10' ] = [ '/WprimeToWZToWhadZinv_narrow_M-3500_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]
	#Samples[ 'QCD_Pt-15to7000' ] = [ '/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 10000 ]


	processingSamples = {}
	if 'all' in args.dataset:
		for sam in Samples: processingSamples[ sam ] = Samples[ sam ]
	else:
		for sam in Samples:
			if sam.startswith( args.dataset ): processingSamples[ sam ] = Samples[ sam ]
	if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'


	for sam in processingSamples:
		dataset = processingSamples[sam][0]

		config.Data.inputDataset = dataset
		config.Data.unitsPerJob = processingSamples[sam][1]
		config.Site.storageSite = args.storageSite
		config.Data.outLFNDirBase = '/store/user/'+os.environ['USER']+'/JetObservables/'+args.version

		if 'BTagMu' in dataset:
			procName = dataset.split('/')[1]+'_'+dataset.split('/')[2]+'_'+args.version
			config.Data.lumiMask = ( 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt' if '2017' in sam else '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' )
			config.Data.splitting = 'LumiBased'
			config.General.workArea = 'crab_projects/Data'

		else:
			procName = dataset.split('/')[1]+('Ext' if 'Ext' in sam else '')+'_'+args.version
			config.Data.splitting = 'EventAwareLumiBased'
			config.General.workArea = 'crab_projects/'

		#config.JobType.pyCfgParams = listParam
		config.General.requestName = procName

		print config
		print '|--- Submmiting sample: ', procName
		if not args.testNoSend: submit(config)
