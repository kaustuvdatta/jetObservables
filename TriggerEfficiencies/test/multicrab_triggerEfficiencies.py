#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def submitJobs( job, inputFiles, unitJobs ):

    from CRABAPI.RawCommand import crabCommand
    from WMCore.Configuration import Configuration
    config = Configuration()

    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.                                                        =
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.section_("General")
    config.General.workArea = options.dir
    #config.General.transferLogs = False
    #config.General.transferOutputs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'triggerEfficiencies_MiniAOD.py'
    config.JobType.allowUndistributedCMSSW = True

    config.section_("Data")
    #config.Data.publication = True
    #config.Data.publishDBS = 'phys03'
    #config.Data.inputDBS = 'phys03'
    #config.Data.ignoreLocality = True

    config.section_("Site")
    config.Site.storageSite = options.storageSite
    ##config.Site.whitelist = ['T1_US_FNAL','T2_CH_CSCS','T3_US_FNALLPC']
    #config.Site.blacklist = ['T2_US_Florida','T3_TW_*','T2_BR_*','T2_GR_Ioannina','T2_BR_SPRACE','T2_RU_IHEP','T2_PL_Swierk','T2_KR_KNU','T3_TW_NTU_HEP']


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers


    requestname = 'jetObservables_triggerEfficiencies_'+ job + '_' +options.version
    print requestname

    if job.startswith(('Single', 'JetHT')): config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    config.Data.inputDataset = inputFiles
    config.Data.splitting = 'Automatic'
    #config.Data.unitsPerJob = unitJobs
    #config.Data.outputPrimaryDataset = job

    # since the input will have no metadata information, output can not be put in DBS
    config.Data.outLFNDirBase = '/store/user/'+os.environ['USER']+'/tmpFiles/jetObservables/triggerEfficiencies'

    if len(requestname) > 100: requestname = (requestname[:95-len(requestname)])
    print 'requestname = ', requestname
    config.General.requestName = requestname
    config.Data.outputDatasetTag = requestname
    print 'Submitting ' + config.General.requestName + ', dataset = ' + job
    print 'Configuration :'
    print config
    submit(config)
    #try : submit(config)
    #except : print 'Not submitted.'



if __name__ == '__main__':

    usage = ('usage: python multicrab_triggerEfficiencies.py --datasets NAMEOFDATASET -d DIR -v VERSION')

    parser = OptionParser(usage=usage)
    parser.add_option(
            "-D", "--dir",
            dest="dir", default="crab_projects",
            help=("The crab directory you want to use "),
            )
    parser.add_option(
            "-d", "--datasets",
            dest="datasets", default='all',
            help=("File listing datasets to run over"),
            )
    parser.add_option(
            "-S", "--storageSite",
            dest="storageSite", default="T3_CH_PSI",
            help=("Storage Site"),
            )
    parser.add_option(
            "-v", "--version",
            dest="version", default="102X_v00",
            help=("Version of output"),
            )
    parser.add_option(
            "-p", "--pythonFile",
            dest="pythonFile", default="triggerEfficiencies.py",
            help=("python file to run"),
            )


    (options, args) = parser.parse_args()


    dictSamples = {}

    dictSamples['JetHT_Run2017B'] = [ '/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD', 1 ]
    dictSamples['JetHT_Run2017C'] = [ '/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD', 1 ]
    dictSamples['JetHT_Run2017D'] = [ '/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD', 1 ]
    dictSamples['JetHT_Run2017E'] = [ '/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD', 1 ]
    dictSamples['JetHT_Run2017F'] = [ '/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD', 1 ]

    processingSamples = {}
    if 'all' in options.datasets:
        for sam in dictSamples: processingSamples[ sam ] = dictSamples[ sam ]
    else:
        for sam in dictSamples:
            if sam.startswith( options.datasets ): processingSamples[ sam ] = dictSamples[ sam ]

    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    for isam in processingSamples:

        if '2017' in isam: options.year = '2017'
        options.datasets = isam

        submitJobs( isam, processingSamples[isam][0], processingSamples[isam][1] )
