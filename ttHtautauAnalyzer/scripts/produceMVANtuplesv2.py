#!/usr/bin/env python
import os
import argparse
from timeit import default_timer as timer
from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getDatasetDict

parser = argparse.ArgumentParser(description='produce MVA ntuples using makeMVANtuple.cc')
parser.add_argument('label', type=str, help="Crab job label for the input event ntuples. E.g. 2019jan")
parser.add_argument('samples', nargs='+', help="sample names")
parser.add_argument('--datasetlist', type=str,
                    default="../dataFiles/DatasetList_2017reMiniAODv2.csv")
parser.add_argument('--analysis', nargs='+',
                    choices=['1l2tau','2lss1tau','3l1tau','2l2tau','control_ttW',
                             'control_ttZ'],
                    default=['1l2tau','2lss1tau','3l1tau','2l2tau','control_ttW',
                             'control_ttZ'],
                    help="analysis type")
parser.add_argument('--jobtypes', nargs='+', choices=['datacard','control'],
                    default=['datacard','control'],
                    help="For making datacards or for control regions")
parser.add_argument('--selection', choices=['application_fake','signal',
                                            'control','control_fakeAR'])
parser.add_argument('-r','--redirector',type=str,
                    default='root://cmsxrootd.fnal.gov/',
                    help="redirector for xrootd")
parser.add_argument('-t','--topeosdir', type=str,
                    default="/store/user/ztao/ttHtaus_94X/")

parser.add_argument('-o','--outdir', type=str, help="Output directory")
#parser.add_argument('-s','--doSystematics', action='store_true',
#                    help="include event weights for systematic uncertainties")
parser.add_argument('-m','--genmatching', action='store_true',
                    help="Require gen matching in event selections")
parser.add_argument('-c','--corrections', nargs='+',
                    choices=['NA','JESUp','JESDown','JERUp','JERDown',
                             'TESUp','TESDown'],
                    default=['NA'], help="Jet/Tau energy correction")
parser.add_argument('--transfer_inputs', action='store_true',
                    help="Copy input files locally from EOS")
parser.add_argument('-l', '--log', type=str, help="Log name")
parser.add_argument('--dryrun', action='store_true',
                    help="Print the command instead of actually executing it")

args = parser.parse_args()

datasets = { 'data_e':'SingleElectron','data_mu':'SingleMuon','data_dieg':'DoubleEG',
             'data_dimu':'DoubleMuon','data_mueg':'MuonEG'}
samplelist_dict = getDatasetDict(args.datasetlist)

logfile = 'mvaNtupleList_'+args.label+'.log'
if args.log is not None:
    logfile = args.log

mvantupleList = open(logfile,'w')

genmatch_str = ' -m true' if args.genmatching else ' -m false'
    
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

    ntuplename += args.label+'/'+filename

    print sample, ntuplename

    # add redirector
    ntuplename = args.redirector + ntuplename
    
    if args.transfer_inputs:
        if args.dryrun:
            print 'xrdcp -s '+ntuplename+' .'
        else:
            os.system('xrdcp -s '+ntuplename+' .')
        ntuplename = filename
        
    for atype in args.analysis:
        print atype
        
        anatype = atype
        # selection types
        seltypes = ['signal_'+anatype, 'control_'+anatype]

        if 'data' in sample:  # collision data
            seltypes.append('application_fake_'+anatype)
            seltypes.append('control_fakeAR_'+anatype)
            if anatype=='2lss1tau':
                seltypes.append("application_flip_2lss1tau")
                seltypes.append("control_flipAR_2lss1tau")
                
        # special cases
        if atype=='control_ttW':
            anatype = '2lss'
            seltypes = ['control_ttW']
            if 'data' in sample:
                seltypes += ['control_fakeAR_ttW','control_flipAR_ttW']
        if atype=='control_ttZ':
            anatype = '3l'
            seltypes = ['control_ttZ']
            if 'data' in sample:
                seltypes.append('control_fakeAR_ttZ')
        if atype=='control_WZ':
            anatype = '3l'
            seltypes = ['control_WZ']
            if 'data' in sample:
                seltypes.append('control_fakeAR_WZ')

        if args.selection is not None: # if selection types are explicitly set
            seltypes = [args.selection+'_'+anatype]

            if atype == 'control_ttW':
                if args.selection in ['application_fake', 'control_fakeAR']:
                    seltypes = ['control_fakeAR_ttW','control_flipAR_ttW']
                else:
                    seltypes = ['control_ttW']

            if atype == 'control_ttZ':
                if args.selection in ['application_fake', 'control_fakeAR']:
                    seltypes = ['control_fakeAR_ttZ']
                else:
                    seltypes = ['control_ttZ']

        for seltype in seltypes:

            if 'control_' in seltype:
                if 'control' not in args.jobtypes and args.selection is None:
                    continue
            else:
                if 'datacard' not in args.jobtypes and args.selection is None:
                    continue
            
            print seltype

            if args.outdir is None:
                args.outdir = args.label
            
            if args.dryrun:
                print 'mkdir -p '+args.outdir
            else:
                os.system('mkdir -p '+args.outdir)
                
            print "make output directory:", args.outdir

            outdirectory = args.outdir
            if args.outdir[-1]!='/':
                outdirectory += '/'
                
            outname = outdirectory+'mvaNtuple_'+sample+'_'+seltype+'.root'
            argument = ' -i '+ntuplename+' -o '+outname+' --anatype '+anatype+' --seltype '+seltype+' -s '+sample
            argument+=genmatch_str
            
            if 'data' in sample:
                if args.dryrun:
                    print 'makeMVANtuple'+argument+' --isdata true'
                else:
                    os.system('makeMVANtuple'+argument+' --isdata true')
                mvantupleList.write(outname+'\n')
            else:
                cross_section = samplelist_dict[sample]['xsection']
                
                if cross_section=='' or cross_section is None:
                    cross_section='-99.'

                for ec in args.corrections:
                    if ec.lower() in ['jesup','jesdown','jerup','jerdown','tesup','tesdown']:
                        outname_ec = outname.replace('.root','_'+ec.lower()+'.root')
                        argument_ec = argument.replace(outname, outname_ec)
                        if args.dryrun:
                            print 'makeMVANtuple'+argument_ec+' --isdata false -x '+cross_section+' --correction '+ec
                        else:
                            os.system('makeMVANtuple'+argument_ec+' --isdata false -x '+cross_section+' --correction '+ec)
                        mvantupleList.write(outname_ec+'\n')
                    else: # nominal
                        if args.dryrun:
                            print 'makeMVANtuple'+argument+' --isdata false -x '+cross_section
                        else:
                            os.system('makeMVANtuple'+argument+' --isdata false -x '+cross_section)
                        mvantupleList.write(outname+'\n')

    stop = timer()
    print "Process", sample, "took", (stop-start)/60., "min"

    if args.transfer_inputs:
        if args.dryrun:
            print 'rm '+filename
        else:
            os.system('rm '+filename)
