#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser
from datasets_dijetSel_RunIISummer20UL_allMC_nomWts import dictSamples, checkDict #datasets_dijetSel_RunIISummer20UL_allMC_nomWts

def make_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def createBash():

    BASH_SCRIPT = '''
#this is not meant to be run locally
#
echo Check if TTY
if [ "`tty`" != "not a tty" ]; then
  echo "YOU SHOULD NOT RUN THIS IN INTERACTIVE, IT DELETES YOUR LOCAL FILES"
else

echo "ENV..................................."
env
echo "VOMS"
voms-proxy-info -all
echo "CMSSW BASE, python path, pwd"
echo $CMSSW_BASE
echo $PYTHON_PATH
echo $PWD
rm -rf $CMSSW_BASE/lib/
rm -rf $CMSSW_BASE/src/
rm -rf $CMSSW_BASE/module/
rm -rf $CMSSW_BASE/python/
mv lib $CMSSW_BASE/lib
mv src $CMSSW_BASE/src
mv python $CMSSW_BASE/python

echo Found Proxy in: $X509_USER_PROXY
'''
    open('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+("_isSigMC" if options.isSigMC else "" )+'.sh', 'w').write(BASH_SCRIPT)
    with open('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+("_isSigMC" if options.isSigMC else "" )+'.sh', 'a') as txtfile:
        cmd = "python "+options.pythonFile+" --sample "+options.datasets+" --selection "+options.selection+" --year "+options.year+" --runEra "+options.runEra+(" --onlyUnc "+options.onlyUnc if options.onlyUnc else "" )+(" --onlyTrees" if options.onlyTrees else "" )+(" --isSigMC" if options.isSigMC else "" )+'\nls\nfi'
        txtfile.write(cmd)


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
    config.JobType.psetName = 'PSet.py'
    config.JobType.maxMemoryMB = 6000 #if (options.onlyUnc and options.onlyUnc.startswith('_je')) else 2500
    #if (options.onlyUnc and options.onlyUnc.startswith('_je')): 
    config.JobType.maxJobRuntimeMin = 1200 if (options.onlyUnc and options.onlyUnc.startswith('_jes')) else 800
    config.JobType.numCores = 4
    config.JobType.allowUndistributedCMSSW = True



    config.section_("Data")
    config.Data.inputDBS = 'phys03'
    config.Data.ignoreLocality = True

    config.section_("Site")
    config.Site.storageSite = options.storageSite
    config.Site.whitelist = ['T2_CH_CSCS','T2_CH_CERN','T1_IT_*','T1_FR_*','T1_DE_*','T2_DE_*','T2_IT_*','T2_FR_*'] #,'T2_DE_*','T1_IT_*','T1_FR_*','T2_IT_*','T2_FR_*']#,'T2_HU_*','T1_ES_*','T2_ES_*','T2_PT_*','T1_US_*']


    def submit(config):
        try:
            crabCommand('submit', config = config)#, dryrun=True)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers


    config.JobType.scriptExe = 'runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+("_isSigMC" if options.isSigMC else "" )+'.sh'
    config.JobType.inputFiles = [ options.pythonFile ,'haddnano.py', 'keep_and_drop_dijet.txt']#, 'keep_and_drop.txt']
    #if (options.onlyUnc and options.onlyUnc.startswith('_jes')): config.JobType.maxJobRuntimeMin = 2000 
    config.JobType.sendPythonFolder  = True
    isDataFlag=False
    if job.startswith(('UL17_Single', 'UL17_JetHT', 'UL18_Single', 'UL18_JetHT', 'UL16_Single', 'UL16_JetHT', 'UL16_preVFP_Single', 'UL16_preVFP_JetHT', 'JetHT', 'SingleMuon')):
        isDataFlag=True
        if options.year.startswith('2016'):
            config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
        if options.year.startswith('2017'): 
            config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        elif options.year.startswith('2018'): 
            config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
    config.Data.inputDataset = inputFiles
    
    config.Data.splitting ='LumiBased' if isDataFlag else 'FileBased'#( (('JetHT' in job) or ('1000to1500' in job) or ('2000' in job) or ('500to700' in job)) and options.onlyUnc.startswith('_je') and not ('JetHT' in job)) else 'FileBased'#'Automatic'
    config.Data.unitsPerJob = 50 if config.Data.splitting=='LumiBased' else unitJobs#unitJobs#'Automatic'
    if not(config.Data.splitting)=='LumiBased': config.Data.unitsPerJob = 1 if ('MLM' in job) else unitJobs
    #config.Data.totalUnits = -1
    #config.Data.splitting = 'FileBased'
    #config.Data.unitsPerJob = 2 if (options.onlyUnc and options.onlyUnc.startswith('_jes')) else unitJobs

    # since the input will have no metadata information, output can not be put in DBS
    config.JobType.outputFiles = [ 'jetObservables_nanoskim.root', 'jetObservables_histograms.root']
    config.Data.outLFNDirBase = '/store/user/'+os.environ['USER']+'/jetObservables/'

    requestname = 'jetObs_Skim_'+ job.replace('_','').replace('-','').replace('PSWeights', '')+options.onlyUnc+'_'+options.year+("_isSigMC" if options.isSigMC else "" )+ '_' +options.version
    print requestname
    
    
    if len(requestname) > 100: requestname = (requestname[:95-len(requestname)])
    

    if os.path.isdir('crab_projects/crab_'+requestname):
        print '|-------> JOB '+requestname+' has already a folder. Please remove it.'
        os.remove('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+("_isSigMC" if options.isSigMC else "" )+'.sh')
        return False

    print 'requestname = ', requestname
    config.General.requestName = requestname
    if isDataFlag: config.Data.outputDatasetTag = 'Run'+inputFiles.split('Run')[1].split('v2pt2')[0]+'_jetObsSkim_'+options.onlyUnc+options.version+('EXT'+job.split('EXT')[1] if 'EXT' in job else '')
    else: config.Data.outputDatasetTag = 'Run'+inputFiles.split('Run')[1].split('-106X')[0]+'_jetObsSkim_'+options.onlyUnc+options.version+('EXT'+job.split('EXT')[1] if 'EXT' in job else '')
    print 'Submitting ' + config.General.requestName + ', dataset = ' + job
    print 'Configuration :', options.year+options.runEra if 'JetHT' in job else ''
    print config
    submit(config)
    #try : submit(config)
    #except : print 'Not submitted.'
    os.remove('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+("_isSigMC" if options.isSigMC else "" )+'.sh')



if __name__ == '__main__':

    usage = ('usage: python multicrab_nSubProducer.py --datasets NAMEOFDATASET -d DIR -v VERSION')

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
            "-y", "--year",
            dest="year", default="2017",
            help=("Version of output"),
            )
    parser.add_option(
            "-s", "--selection",
            dest="selection", default="dijet",
            help=("Selection: dijet, Wtop"),
            )
    parser.add_option(
            "-p", "--pythonFile",
            dest="pythonFile", default="jetObservables_nSubProducer_dijetSel.py",
            help=("python file to run"),
            )
    parser.add_option(
            '--runEra',
            action="store",
            help="Run era for data",
            default=""
    )
    parser.add_option(
        '--onlyUnc',
        action="store",
        default='',
        help="Run only specific uncertainty variations"
    )
    parser.add_option(
        '--onlyTrees',
        action="store_true",
        help="Do not save histograms, only trees"
    )
    parser.add_option(
        '--isSigMC',
        action="store_true",
        #default="",
        help="Save branches with sys+puWeights for signal/main MC else store basic reco/gen branches"
    )
   

    (options, args) = parser.parse_args()
    print(options,args)

    processingSamples = {}
    for sam in dictSamples:
        if sam.startswith( options.datasets ) | options.datasets.startswith('all'):
            if checkDict( sam, dictSamples )['selection'] != options.selection: continue
            if sam.startswith(('JetHT', 'SingleMuon')):
                if options.runEra=="":
                    for iera in checkDict( sam, dictSamples )[options.year]['nanoAOD']:
                        #print(sam,iera)
                        processingSamples[ sam+'Run'+options.year+iera ] = [ checkDict( sam, dictSamples )[options.year]['nanoAOD'][iera], 1 ]
                else:
                    for iera in [options.runEra]:
                        #print(sam,iera)
                        processingSamples[ sam+'Run'+options.year+iera ] = [ checkDict( sam, dictSamples )[options.year]['nanoAOD'][iera], 1 ]
                #print(sam,options,args)
            else:
                tmpList = checkDict( sam, dictSamples )[options.year]['nanoAOD']
                processingSamples[ sam ] = [ tmpList[0], 1 ]
                if len(tmpList)>1:
                    for iext in range(1,len(tmpList)):
                        processingSamples[ sam+'EXT'+str(iext) ] = [ tmpList[iext], 1 ]
                options.runEra = ''
        
    #print(processingSamples)
 
    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'


    for isam in processingSamples:
        if 'JetHT' in isam: options.runEra = isam.split(options.year)[1]#iera
        print (isam,options.runEra)
        if not processingSamples[isam][0]:
            print(' Sample ',isam,' does not have nanoAOD stored in datasets.py. Continuing with the next')
            continue
        #if isam.startswith('QCD') or isam.startswith('JetHT'): options.selection = 'dijet'
        #else: options.selection = 'Wtop'

        options.datasets = isam
        print('Creating bash file for %s...'%(isam))
        createBash()

        print ("Submitting dataset %s" % (processingSamples[isam]))#, len(processingSamples[isam][0])))
        submitJobs( isam, processingSamples[isam][0], processingSamples[isam][1] )

