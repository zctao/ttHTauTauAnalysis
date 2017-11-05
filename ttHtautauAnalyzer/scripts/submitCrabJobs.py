#!/usr/bin/env python

# submit multiple crab jobs
# modified from Jorge's script

import os
import argparse
from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import channel_suffix

parser = argparse.ArgumentParser(description='Submit multiple crab jobs for ntuple making')

parser.add_argument('samples', type=str,
                    help="List of MC/data samples to run crab jobs")
parser.add_argument('analysis', choices=['1l2tau','2lss1tau','3l1tau'],
                    help="Analysis type")
parser.add_argument('selection', choices=['signal','loose','fakes','flips'],
                    help="Selection type")

chgroup = parser.add_mutually_exclusive_group(required=True)
chgroup.add_argument('--channels',nargs='+', type=str,
                    help="List of channels to run")
chgroup.add_argument('--channel_list', type=str,
                     help="Text file that stores a list of channels to run")

parser.add_argument('--prefix', type=str, default='',
                    help="Prefix of crab job name.")

parser.add_argument('-c','--config', type=str, default='analyzer2016_cfg.py',
                    help="cmsRun config file")
parser.add_argument('--dbs', choices=['global','phys01','phys02','phys03'],
                    default='global',
                    help="Alias of the URL of the DBS reader instance where the input dataset is published ")
parser.add_argument('--systematics', nargs='+',
                    choices=['jesup','jesdown','tesup','tesdown'],
                    default=[],
                    help="Systematics")
parser.add_argument('--syst_only', action='store_true',
                    help="Only create crab jobs for systematics")
parser.add_argument('-o','--outdir', type=str,
                    default='/store/user/ztao/ttHtautau_80X',
                    help="Output directory")
parser.add_argument('-d','--dryrun', action='store_true',
                    help="Run the script without actually submitting crab jobs")
parser.add_argument('--lumimask', type=str,
                    default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
                    help='Lumi mask to run on data sample')

args = parser.parse_args()

# channels
def getChannelList(infile):
    f=open(infile, 'r')
    tmp=f.read()
    return tmp.replace('\n','').replace(' ','').split(',')
    
channels = args.channels if args.channels is not None else getChannelList(args.channel_list)
#print channels

# sample list
samplelist = args.samples

def getSelectionRegionString(analysis, selection):
    # correspond to the 'SelectionRegion' argument in cmsRun config file
    if analysis == '2lss1tau' and selection=='signal':
        return 'signal_2lss1tau'
    elif analysis == '2lss1tau' and selection=='loose':
        return 'loose_2lss1tau'
    elif analysis == '2lss1tau' and selection=='fakes':
        return 'control_fake_2lss1tau'
    elif analysis == '2lss1tau' and selection=='flips':
        return 'control_2los1tau'
    elif analysis == '1l2tau' and selection=='signal':
        return 'signal_1l2tau'
    elif analysis == '1l2tau' and selection=='loose':
        return 'loose_1l2tau'
    elif analysis == '1l2tau' and selection=='fakes':
        return 'control_fake_1l2tau'
    elif analysis == '3l1tau' and selection=='signal':
        return 'signal_3l1tau'
    #elif analysis_tpye == '3l1tau' and selection=='loose':
    #    return 'loose_3l1tau'
    elif analysis == '3l1tau' and selection=='fakes':
        return 'control_fake_3l1au'
    else:
        assert(False)
        return ''

    
# CRAB config file template
template = '''from CRABClient.UserUtilities import config
config = config()
config.General.requestName = '%(name)s'
config.General.transferLogs = True
config.General.workArea = '/uscms/home/ztao/nobackup/crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(configfile)s'
config.JobType.pyCfgParams = %(cfgparams)s
config.JobType.sendExternalFolder = True
config.Data.inputDataset = %(dataset)s
config.Data.inputDBS = '%(DBS)s'
config.Data.splitting = '%(splitting)s'
config.Data.unitsPerJob = %(unit)s
%(lumimask)s
%(runrange)s
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '%(outdir)s'
config.Site.storageSite = 'T3_US_FNALLPC'
'''

#'analyzer2016_cfg.py'

for ch in channels:

    print ch
    
    dataset = ''
    pset = []
    jobnames = []
    
    with open(samplelist) as f:
        for line in f:

            if line.strip()!=ch:
                    continue
            
            if not 'data' in ch: # MC samples

                dataset = f.next().strip()
                perjob = f.next().strip()
                    
                ## determine parameter set for the job
                pset0 = "['SampleName=','AnalysisType=','SelectionRegion=','doCutFlow=True', 'doSystematics=True','JECType=NA','TauESType=NA']"
                    
                # sample name
                sname = ch
                for suffix in channel_suffix:
                    sname=sname.replace(suffix,'')
                pset0=pset0.replace("SampleName=","SampleName="+sname)
                    
                # analysis type
                pset0=pset0.replace("AnalysisType=","AnalysisType="+args.analysis)
                # selection type
                pset0=pset0.replace("SelectionRegion=","SelectionRegion="+getSelectionRegionString(args.analysis, args.selection))
                
                if not args.syst_only:
                    pset.append(pset0)
                    jobnames.append(ch+'_'+args.analysis)
                    
                # systematicss
                for syst in args.systematics:
                    # turn off other systematics computations not related to jet or tau energy scale
                    pset_syst = pset0.replace("doSystematics=True","doSystematics=False")
                    if syst=='jesup':
                        pset_syst = pset_syst.replace("JECType=NA","JECType=JESUp")
                    elif syst=='jesdown':
                        pset_syst = pset_syst.replace("JECType=NA","JECType=JESDown")
                    elif syst=='tesup':
                        pset_syst = pset_syst.replace("TauESType=NA","TauESType=tauESUp")
                    elif syst=='tesdown':
                        pset_syst = pset_syst.replace("TauESType=NA","TauESType=tauESDown")

                    pset.append(pset_syst)
                    jobnames.append(ch+'_'+syst+'_'+args.analysis)

                assert(len(pset)==len(args.systematics) if args.syst_only else len(args.systematics)+1)
                assert(len(jobnames)==len(args.systematics) if args.syst_only else len(args.systematics)+1)
                
            else:  # collision data

                dataset = f.next().strip()
                run = f.next().strip()
                perjob = f.next().strip()

                pset0="['isData=True','SampleName=','AnalysisType=','SelectionRegion=','TurnOffHLTCut=True', 'HIPSafeMediumMuon=True', 'Is2016H=False']"

                # sample name
                sname = ''
                if 'signal' in args.selection:
                    sname = 'data_obs'
                elif 'fakes' in args.selection:
                    sname = 'fakes_data'
                elif 'flips' in args.selection:
                    sname = 'flips_data'
                else:
                    print 'Selection argument: ', args.selection, ' is not supported for data sample'
                    exit()

                pset0=pset0.replace("SampleName=","SampleName="+sname)

                # analysis type
                pset0=pset0.replace("AnalysisType=","AnalysisType="+args.analysis)
                # selection type
                pset0=pset0.replace("SelectionRegion=","SelectionRegion="+getSelectionRegionString(args.analysis, args.selection))

                if '_2016g' in ch or '_2016h' in ch:
                    pset0=pset0.replace("HIPSafeMediumMuon=True","HIPSafeMediumMuon=False")
                if '_2016h' in ch:
                    pset0=pset0.replace("Is2016H=False","Is2016H=True")

                pset.append(pset0)
                jobnames.append(ch.replace('data',sname)+'_'+args.analysis)
                assert(len(pset)==1)
                assert(len(jobnames)==1)

    # write crab config files for each of the pset
    for ps, jn in zip(pset, jobnames):
        
        vd = locals()
        vd['name'] = args.prefix+jn
        vd['dataset'] = dataset
        vd['configfile'] = args.config
        vd['cfgparams'] = ps
        vd['DBS'] = args.dbs
        vd['unit'] = perjob
        vd['outdir'] = args.outdir
        if 'data' in ch:
            vd['lumimask'] = "config.Data.lumiMask = "+"'"+args.lumimask+"'"
            vd['runrange'] = 'config.Data.runRange = '+ run
            vd['splitting'] = 'LumiBased'
        else:
            vd['lumimask'] = ''
            vd['runrange'] = ''
            vd['splitting'] = 'EventAwareLumiBased'

        # write crab config files
        open('crab/crabConfig_'+vd['name']+'.py','wt').write(template % vd)
        
        # submit crab jobs
        if not args.dryrun:
            os.system('crab submit -c crab/crabConfig_'+vd['name']+'.py')
