#####################################
#####   This script is meant to run in t3ui02 (PSI)
#####   Requirements:
#####    * uproot: download miniconda and then
#####      conda config --add channels conda-forge   # if you haven't added conda-forge already
#####      conda install uproot
#####################################


#!/usr/bin/env python
import argparse, os, shutil
import numpy as np
import uproot

def convertRootToNumpy( inputFile, outputFolder, outputName ):
    """docstring for convertRootToNumpy"""

    fullTree = uproot.open(inputFile)["addingNsub/SimpleGentree"]
    arrays = fullTree.arrays([
            "genJetPt",
            "genJetMass",
            "genJetEta",
            "genJetPhi",
            "genJetEnergy",
            "genJetNsub0p5tau1",
            "genJetNsub0p5tau2",
            "genJetNsub0p5tau3",
            "genJetNsub0p5tau4",
            "genJetNsub0p5tau5",
            "genJetNsub0p5tau6",
            "genJetNsub0p5tau7",
            "genJetNsub1tau1",
            "genJetNsub1tau2",
            "genJetNsub1tau3",
            "genJetNsub1tau4",
            "genJetNsub1tau5",
            "genJetNsub1tau6",
            "genJetNsub1tau7",
            "genJetNsub2tau1",
            "genJetNsub2tau2",
            "genJetNsub2tau3",
            "genJetNsub2tau4",
            "genJetNsub2tau5",
            "genJetNsub2tau6",
            "genJetNsub2tau7",
            ])

    ## Convert the tree to a npy
    np.save(outputFolder+'/'+outputName, arrays)
    print('Done')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--inputFolder", action='store', dest="inputFolder", default="", help="Path of the folder where the root files are located (wihtout the root://blah//)." )
    parser.add_argument("-o", "--outputFolder", action='store', dest="outputFolder", default="test", help="Name of output folder" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dictFolders = {}
    if args.inputFolder:
        dictFolders[ args.outputFolder ] = [ args.inputFolder ]
    else:
        args.inputFolder = "/pnfs/psi.ch/cms/trivcat/store/user/algomez/jetObservables_"+args.version+"/CRAB_UserFiles/"
        args.outputFolder = args.outputFolder+'/'+args.version+'/'
        dictFolders[args.outputFolder+"QCDPt170to300"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_']
        dictFolders[args.outputFolder+"QCDPt300to470"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_300to470' ]
        dictFolders[args.outputFolder+"QCDPt470to600"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_470to600' ]
        dictFolders[args.outputFolder+"QCDPt600to800"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_600to800' ]
        dictFolders[args.outputFolder+"QCDPt800to1000"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_800to1000' ]
        dictFolders[args.outputFolder+"QCDPt1000to1400"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_1000to1400' ]
        dictFolders[args.outputFolder+"QCDPt1400to1800"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_1400to1800' ]
        dictFolders[args.outputFolder+"QCDPt1800to2400"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_1800to2400' ]
        dictFolders[args.outputFolder+"QCDPt2400to3200"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_2400to3200' ]
        dictFolders[args.outputFolder+"QCDPt3200toInf"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_QCD_Pt_3200toInf' ]
        dictFolders[args.outputFolder+"TTJets"] = [ args.inputFolder+'/jetObservablesInput_v00_94X_TTJets' ]

    for isample, ifolder  in dictFolders.items():

        if not os.path.exists( isample ):
            os.makedirs(isample)
            print("Creating ", isample, " folder.")
        else:
            print("Removing and re-creating ", isample, " folder.")
            shutil.rmtree(isample)
            os.makedirs(isample)

        for root, dirs, allfiles in os.walk( ifolder[0] ):
            for ifile in allfiles:
                if not '.root' in ifile: continue
                if 'histos' in ifile: continue
                print('Converting ', root+"/"+ifile)
                convertRootToNumpy( "root://t3dcachedb03.psi.ch/"+root+'/'+ifile, isample, ifile.split(".root")[0] )

        print("Done. Output files location: ", isample)
