#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser
from datasets_WtopSel_RunIISummer20UL import dictSamples, checkDict

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
###ls -lR .
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
#ls
#echo "python {pythonFile} --sample {datasets} --selection {selection} --year {year}"
#python {pythonFile} --sample {datasets} --selection {selection} --year {year} --runEra {runEra} --onlyTrees
#python {pythonFile} --sample {datasets} --selection {selection} --year {year} --runEra {runEra} 
#python {pythonFile} --sample {datasets} --selection {selection} --year {year} --onlyUnc {onlyUnc} --runEra {runEra}  
#fi
    '''
    #open('runPostProc'+options.datasets+options.year+options.onlyUnc+options.selection+'.sh', 'w').write(BASH_SCRIPT.format(**options.__dict__))
    open('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+options.selection+'.sh', 'w').write(BASH_SCRIPT)
    with open('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+options.selection+'.sh', 'a') as txtfile:
        cmd = "python "+options.pythonFile+" --sample "+options.datasets+" --selection "+options.selection+" --year "+options.year+" --runEra "+options.runEra+(" --onlyUnc "+options.onlyUnc if options.onlyUnc else "" )+(" --onlyTrees" if options.onlyTrees else "" )+(" --isSigMC" if options.isSigMC else "" )+(" --isBkgMC" if options.isBkgMC else "" )+'\nls\nfi'
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
    config.JobType.maxJobRuntimeMin = 600
    config.JobType.numCores = 4
    config.JobType.allowUndistributedCMSSW = True


    config.section_("Data")
    config.Data.inputDBS = 'phys03'
    config.Data.ignoreLocality = True

    config.section_("Site")
    config.Site.storageSite = options.storageSite
    config.Site.whitelist = ['T2_CH_CSCS','T2_CH_CERN','T1_IT_*','T1_FR_*','T1_DE_*','T2_DE_*','T2_IT_*','T2_FR_*']


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers


    config.JobType.scriptExe = 'runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+options.selection+'.sh'
    config.JobType.inputFiles = [ options.pythonFile ,'haddnano.py', 'keep_and_drop_Wtop.txt']
    config.JobType.sendPythonFolder  = True
    
    isMC = True
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

    #config.Data.splitting = 'EventAwareLumiBased' if job.startswith('QCD_Pt') else 'FileBased'
    #config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'EventAwareLumiBased' #if job.startswith('QCD_HT') else 'FileBased'
    #config.Data.unitsPerJob = 900
    config.Data.splitting = 'FileBased' if not('SingleMuon' in job) else 'LumiBased'
    config.Data.unitsPerJob = unitJobs if config.Data.splitting=='FileBased' else 50

    # since the input will have no metadata information, output can not be put in DBS
    config.JobType.outputFiles = [ 'jetObservables_nanoskim.root', 'jetObservables_histograms.root']
    config.Data.outLFNDirBase = '/store/user/'+os.environ['USER']+'/jetObservables/'

    requestname = 'jetObs_Skimmer_'+ job.replace('_','').replace('-','')+options.onlyUnc+'_'+options.year +'_'+options.selection + '_' +options.version

    if len(requestname) > 100: requestname = (requestname[:-(len(requestname)-100)])

    if os.path.isdir('crab_projects/crab_'+requestname):
        print '|-------> JOB '+requestname+' has already a folder. Please remove it.'
        os.remove('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+options.selection+'.sh')
        return False

    print 'requestname = ', requestname
    isvNFlag=True if 'v2pt' in inputFiles else False
    config.General.requestName = requestname

    if (isDataFlag and isvNFlag) or (isvNFlag and 'SingleMuon' not in job) :#('QCD' in job or 'MuEn'in job) and (options.selection.startswith('_W') or options.selection.startswith('_top')):
        s=['v2pt2','v2pt3','v2pt4']
        flag=""
        for x in s:
            if x in inputFiles: 
                flag=x
        #flag = next((x for x in s if x in inputFiles), "")#x for x in s if x in inputFiles else ""
        #print (s,flag)
        if not(flag==""): 
            if 'APV' in inputFiles and options.year.endswith('VFP') :
                flag='APV'+flag
            config.Data.outputDatasetTag = 'Run'+inputFiles.split('Run')[1].split(flag)[0]+flag+'_jetObsSkim_'+options.onlyUnc+options.version+('EXT'+job.split('EXT')[1] if 'EXT' in job else '')
        else:
            k='PFNanoAOD-' if 'PFNanoAOD-' in inputFiles else 'PFNano-'
            config.Data.outputDatasetTag = 'Run'+inputFiles.split('Run')[1].split(k)[0]+'jetObsSkim_'+options.onlyUnc+options.version+('EXT'+job.split('EXT')[1] if 'EXT' in job else '')#.outputDatasetTag
    else: 
        k='PFNanoAOD-' if 'PFNanoAOD-' in inputFiles else 'PFNano-'
        config.Data.outputDatasetTag = 'Run'+inputFiles.split('Run')[1].split(k)[0]+'jetObsSkim_'+options.onlyUnc+options.version+('EXT'+job.split('EXT')[1] if 'EXT' in job else '')#.outputDatasetTag

    print 'Submitting ' + config.General.requestName + ', dataset = ' + job
    print ('inputFiles=',inputFiles, 'outputDatasetTag=', config.Data.outputDatasetTag)
    print 'Configuration :', options.year+options.runEra if 'SingleMuon' in job else ''

    print config
    submit(config)
    #try : submit(config)
    #except : print 'Not submitted.'
    os.remove('runPostProc'+options.datasets+options.onlyUnc+'_'+options.year+options.runEra+options.selection+'.sh')



if __name__ == '__main__':

    usage = ('usage: python multicrab_nSubProducer.py --datasets NAMEOFDATASET -d DIR -v VERSION -y YEAR --onlyUnc _UNC ')

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
            dest="year", default='2017',
            help=("Sample year"),
            )
    parser.add_option(
            "-s", "--selection",
            dest="selection", default="_WtopSel",
            help=("Selection: _WtopSel, _WSel, _topSel"),
            )
    parser.add_option(
            "-p", "--pythonFile",
            dest="pythonFile", default="jetObservables_nSubProducer_WtopSel.py",
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
        help="Save branches with sys+puWeights for signal/main MC else store basic reco/gen branches"
            )  
    parser.add_option(
        '--isBkgMC',
        action="store_true",
        help="Save branches with basic reco/gen branches for bkg MC and signal variation models"
            )  

    (options, args) = parser.parse_args()
    print (options.selection)
    if options.selection.startswith(('_Wtop','Wtop')): sel = 'Wtop'
    elif options.selection.startswith(('_WSel','WSel')): sel = 'W'
    elif options.selection.startswith(('_top','top','_topSel')): sel = 'top'
    processingSamples = {}
    for sam in dictSamples:
        if sam.startswith( options.datasets ) | options.datasets.startswith('all'):
            if checkDict( sam, dictSamples )['selection'] != sel:   
                print("wrong sel",sel)  
                continue
            if sam.startswith(('SingleMuon')):
                if options.runEra=="":
                    for iera in checkDict( sam, dictSamples )[options.year]['nanoAOD']:
                        processingSamples[ sam+'Run'+options.year+iera ] = [ checkDict( sam, dictSamples )[options.year]['nanoAOD'][iera], 1 ]
                        #options.runEra = iera
                        #print (iera,options.runEra)
                else:
                    for iera in [options.runEra]:
                        processingSamples[ sam+'Run'+options.year+iera ] = [ checkDict( sam, dictSamples )[options.year]['nanoAOD'][iera], 1 ]
            else:
                tmpList = checkDict( sam, dictSamples )[options.year]['nanoAOD']
                processingSamples[ sam ] = [ tmpList[0], 1 ]
                if len(tmpList)>1:
                    for iext in range(1,len(tmpList)):
                        processingSamples[ sam+'EXT'+str(iext) ] = [ tmpList[iext], 1 ]
                options.runEra=""

    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    for isam in processingSamples:
        if 'SingleMuon' in isam: options.runEra=isam.split(options.year)[1]
        print(isam,options.runEra)
        if not processingSamples[isam][0]:
            print(' Sample ',isam,' does not have nanoAOD stored in datasets.py. Continuing with the next')
            continue
        
        options.datasets = isam
        print('Creating bash file for %s...'%(isam))
        createBash()

        print ("Submitting dataset %s" % (processingSamples[isam]))
        submitJobs( isam, processingSamples[isam][0], processingSamples[isam][1] )
