#!/usr/bin/env python

import os
import argparse
import ttHTauTauAnalysis.ttHtautauAnalyzer.datacards as dc
from ttHTauTauAnalysis.ttHtautauAnalyzer.cross_section import CrossSection

parser = argparse.ArgumentParser(description='produce MVA ntuples using makeMVANtuple.cc')

parser.add_argument('anatype', choices=['1l2tau','2lss1tau','3l1tau'],
                    help='analysis type')
parser.add_argument('samples', nargs='+',
                    choices=['ttH','TTW','TTZ','TTGJets','TGJets','WG','ZG','WZ','ZZ','WW','WWds','WpWp','WZZ','WWZ','WWW','ZZZ','tZq','TTTT','fakes_data','flips_data','data_obs'],
                    help="sample names")
parser.add_argument('-d','--datasets', nargs='+', choices=['SingleMuon','SingleElectron','DoubleMuon','DoubleEG','MuonEG'],
                    default=['SingleMuon','SingleElectron','DoubleMuon','DoubleEG','MuonEG'], help="dataset for collision data sample")
parser.add_argument('-c','--corrections', nargs='+',
                    choices=['jesup','jesdown','tesup','tesdown'], default=None,
                    help="Jet/Tau energy correction")
parser.add_argument('-s','--doSystematics', action='store_true',
                    help="include event weights for systematic uncertainties")
parser.add_argument('-o','--outdir', type=str, default='./', help="Output directory")
parser.add_argument('-l','--list', type=str, default='ntuplelist.log',
                    help="list of event ntuple to be processed")
parser.add_argument('-r','--redirector', type=str,
                    default='root://cmsxrootd.fnal.gov/',help="redirector for xrootd")
parser.add_argument('-v','--version', choices=['2016','2017','test'],
                    help="Version of analysis")

args = parser.parse_args()

mvantupleList = open(args.outdir+'mvaNtuples.txt','w')

dosys_str = 'true' if args.doSystematics else 'false'

for sample in args.samples:

    # determine seltype
    seltype = 'signal_'+args.anatype
    if 'fakes_data' in sample:
        seltype = 'control_fake_'+args.anatype
    elif 'flips_data' in sample:
        assert(anatype=='2lss1tau')
        seltype = 'control_2los1tau'
    
    if 'data' in sample:
        assert(len(args.datasets)>0)
        for dataset in args.datasets:
            fname = dc.getNtupleFileName_data(args.list, args.anatype, sample, dataset)
            if fname is None:
                continue
            
            ntuplename = args.redirector+fname
            outputname = args.outdir+'mvaNtuple_'+sample+'_'+dataset+'_'+args.anatype+'.root'
            
            os.system('makeMVANtuple -i '+ntuplename+' -o '+outputname+' --anatype '+args.anatype+' --seltype '+seltype+' -d true -m false -u true -w true -s '+dosys_str+' -v '+args.version)
            mvantupleList.write(outputname+'\n')
    else:
        fname = dc.getNtupleFileName_mc(args.list, args.anatype, sample)
        if fname is None:
            continue

        xsection = CrossSection[sample]
        
        ntuplename = args.redirector+fname
        outputname = args.outdir+'mvaNtuple_'+sample+'_'+args.anatype+'.root'
        os.system('makeMVANtuple -i '+ntuplename+' -o '+outputname+' --anatype '+args.anatype+' --seltype '+seltype+' -u true -w true -s '+dosys_str+' -v '+args.version+' -x '+str(xsection))
        mvantupleList.write(outputname+'\n')

        if args.doSystematics:
            for cor in args.corrections:
                fname_cor = dc.getNtupleFileName_mc(args.list, args.anatype, sample, correction=cor)
                if fname_cor is None:
                    continue
                
                ntuplename = args.redirector+fname_cor
                outputname = args.outdir+'mvaNtuple_'+sample+'_'+cor+'_'+args.anatype+'.root'
                
                os.system('makeMVANtuple -i '+ntuplename+' -o '+outputname+' --anatype '+args.anatype+' --seltype '+seltype+' -u true -w true -s '+dosys_str+' -v '+args.version+' -x '+str(xsection))
                mvantupleList.write(outputname+'\n')
