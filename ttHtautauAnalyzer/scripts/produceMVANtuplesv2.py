#!/usr/bin/env python
import os
import argparse
from timeit import default_timer as timer
from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getDatasetDict

parser = argparse.ArgumentParser(description='produce MVA ntuples using makeMVANtuple.cc')

parser.add_argument('samples', nargs='+', help="sample names")
parser.add_argument('--datasetlist', type=str,
                    default="../dataFiles/DatasetList_2017reMiniAODv2.csv")
parser.add_argument('--analysis', nargs='+',
                    choices=['1l2tau','2lss1tau','3l1tau','2l2tau'],
                    help="analysis type")
parser.add_argument('--selection', nargs='+', type=str, help="selection type")
parser.add_argument('-r','--redirector',type=str,
                    default='root://cmsxrootd.fnal.gov/',
                    help="redirector for xrootd")
parser.add_argument('-t','--topeosdir', type=str,
                    default="/store/user/ztao/ttHtaus_94X/")
parser.add_argument('--version', type=str, default="may2018",
                    help="Event ntuple version")
parser.add_argument('-o','--outdir', type=str, default="./", help="Output directory")
#parser.add_argument('-s','--doSystematics', action='store_true',
#                    help="include event weights for systematic uncertainties")
parser.add_argument('-c','--corrections', nargs='+',
                    choices=['jesup','jesdown','tesup','tesdown'], default=[],
                    help="Jet/Tau energy correction")
parser.add_argument('--transfer_inputs', action='store_true',
                    help="Copy input files locally from EOS")

args = parser.parse_args()

datasets = { 'data_e':'SingleElectron','data_mu':'SingleMuon','data_dieg':'DoubleEG',
             'data_dimu':'DoubleMuon','data_mueg':'MuonEG'}
samplelist_dict = getDatasetDict(args.datasetlist)

mvantupleList = open('mvaNtupleList_'+args.version+'.log','w')

anatypes=['1l2tau','2lss1tau','3l1tau','2l2tau']
if args.analysis is not None:  # if analysis types are explicitly set
    anatypes = args.analysis
    
for sample in args.samples:
    start = timer()
    
    # get event ntuple name
    ntuplename = args.topeosdir
    if 'data' in sample:
        assert(sample in datasets)
        ntuplename += datasets[sample]+'/'
    else:
        assert(sample in samplelist_dict)
        ntuplename += samplelist_dict[sample]['dataset'].split('/')[1]+'/'

    filename = 'ntuple_'
    if 'data' in sample:
        filename += 'data_incl.root'
    else:
        filename += sample+'_incl.root'

    ntuplename += args.version+'/'+filename

    print sample, ntuplename

    # add redirector
    ntuplename = args.redirector + ntuplename
    
    if args.transfer_inputs:
        os.system('xrdcp -s '+ntuplename+' .')
        ntuplename = filename
        
    for anatype in anatypes:
        print anatype
        
        # selection types
        seltypes = ['signal_'+anatype, 'control_fake_'+anatype]
        if anatype=='2lss1tau':
            seltypes.append('control_ttW')
        if anatype=='3l1tau':
            seltypes.append('control_ttZ')
            #seltypes.append('control_WZ')  # TODO
            
        if 'data' in sample:  # collision data
            seltypes.append('application_fake_'+anatype)
            seltypes.append('control_fakeAR_'+anatype)
            if anatype=='2lss1tau':
                seltypes.append("application_flip_2lss1tau")

        if args.selection is not None: # if selection types are explicitly set
            seltypes = args.selection

        for seltype in seltypes:
            print seltype

            os.system('mkdir -p '+args.outdir+args.version+'/')
            print "make output directory:", args.outdir+args.version+'/'
            
            outname = args.outdir+args.version+'/mvaNtuple_'+sample+'_'+seltype+'.root'
            argument = ' -i '+ntuplename+' -o '+outname+' --anatype '+anatype+' --seltype '+seltype
            
            if 'data' in sample:
                os.system('makeMVANtuple'+argument+' --isdata true')
                mvantupleList.write(outname+'\n')
            else:
                cross_section = samplelist_dict[sample]['xsection']
                os.system('makeMVANtuple'+argument+' --isdata false -x '+cross_section)
                mvantupleList.write(outname+'\n')
                
                for ec in args.corrections:
                    ntuplename_ec = ntuplename.replace('_incl.root',
                                                       '_incl_'+ec+'.root')
                    outname_ec = outname.replace('.root','_'+ec+'.root')
                    argument_ec = argument.replace(ntuplename, ntuplename_ec)
                    argument_ec = argument.replace(outname, outname_ec)
                    
                    os.system('makeMVANtuple'+argument+' --isdata false -x '+cross_section)
                    mvantupleList.write(outname_ec+'\n')

    stop = timer()
    print "Process", sample, "took", (stop-start)/60., "min"

    if args.transfer_inputs:
        os.system('rm '+filename)
