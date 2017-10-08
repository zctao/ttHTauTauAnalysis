# submit multiple crab jobs
# modified from Jorge's script

import os, sys

import argparse
parser = argparse.ArgumentParser(description='Submit multiple crab jobs for ntuple making')

parser.add_argument('-s', '--samples', type=str,
                    default='../data/SampleList_FastSim2016.txt',
                    help="List of MC/data samples to run crab jobs")
parser.add_argument('-c', '--channels', choices=['train','2016','2017'],
                    default='train',
                    help="Lists of job names: saved in crab_channels.py. Edit this module if need to add/change channels to run.")
parser.add_argument('-a','--analysis', choices=['1l2tau','2lss1tau','3l1tau'],
                    default='2lss1tau', help="Analysis type")
parser.add_argument('-l','--loose', action='store_true',
                    help='loose selection')
parser.add_argument('-o','--outdir', type=str,
                    default='/store/user/ztao/ttH_80X',
                    help="Output directory")
parser.add_argument('--selection', type=str, default='',
                    help="Selection type is normally determined by the sample name and analysis type. This argument, if not empty string, overwrites the usual selection type in case of special need.")
parser.add_argument('-d','--dryrun', action='store_true',
                    help='Run the script without actually submitting crab jobs')

args = parser.parse_args()

import crab_channels
channels = []
if args.channels=='train':
    channels = crab_channels.channels_train
elif args.channels=='2016':
    channels = crab_channels.channels_2016
elif args.channels=='2017':
    channels = crab_channels.channels_2017

samplelist = args.samples

string = '''from CRABClient.UserUtilities import config
config = config()
config.General.requestName = '%(name)s'
config.General.transferLogs = True
config.General.workArea = '/uscms/home/ztao/nobackup/crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyzer2016_cfg.py'
config.JobType.pyCfgParams = %(cfgparams)s
config.JobType.sendExternalFolder = True
config.Data.inputDataset = %(dataset)s
config.Data.inputDBS = 'phys03'
config.Data.splitting = '%(splitting)s'
config.Data.unitsPerJob = %(unit)s
%(lumimask)s
%(runrange)s
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '%(outdir)s'
config.Site.storageSite = 'T3_US_FNALLPC'
'''

samplename_suffix = ['_ext','_jesup','_jesdown','_tesup','_tesdown',
                     '_p1','_p2','_p3']

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]

sample=""
pset=""

for ch in channels:
    
    del sample
    del pset

    print ch
    
    with open(samplelist) as f:
        for line in f:
            if not 'data' in ch:
                if line.strip()==ch.replace("_jesup","") or line.strip()==ch.replace("_jesdown","") or line.strip()==ch.replace("_tesup","") or line.strip()==ch.replace("_tesdown",""):
                    sample = f.next().strip()
                    perjob = f.next().strip()
                    
                    ## determine parameter set for the job
                    pset = "['SampleName=','AnalysisType=','SelectionRegion=','doCutFlow=True', 'doSystematics=True','JECType=NA','TauESType=NA']"
                    
                    # sample name
                    sname = ch
                    for suffix in samplename_suffix:
                        sname=sname.replace(suffix,'')
                    pset=pset.replace("SampleName=","SampleName="+sname)
                    # analysis type
                    pset=pset.replace("AnalysisType=","AnalysisType="+args.analysis)
                    # selection type
                    stype='loose_'+args.analysis if args.loose else 'signal_'+args.analysis
                    if args.selection != '':
                        stype=args.selection
                    pset=pset.replace("SelectionRegion=","SelectionRegion="+stype)

                    # systematicss
                    if '_jesup' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESUp")
                    if '_jesdown' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("JECType=NA","JECType=JESDown")

                    if '_tesup' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("TauESType=NA","TauESType=tauESUp")

                    if '_tesdown' in ch:
                        pset=pset.replace("doSystematics=True","doSystematics=False")
                        pset=pset.replace("TauESType=NA","TauESType=tauESDown")
                        
                    
                    break
                
            else:
                if line.strip()==remove_prefix(ch,"data_obs") or line.strip()==remove_prefix(ch,"fakes_data") or line.strip()==remove_prefix(ch,"flips_data"):
                    sample = f.next().strip()
                    run = f.next().strip()
                    perjob = f.next().strip()

                    ## FIXME for 1l2tau and 3l1tau 
                    
                    pset="['isData=True','SampleName=','SelectionRegion=','TurnOffHLTCut=True', 'HIPSafeMediumMuon=True', 'Is2016H=False']"
                    if 'fakes_' in ch:
                        pset=pset.replace("SampleName=","SampleName=fakes_data")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=control_fake_2lss1tau")
                    elif 'flips_' in ch:
                        pset=pset.replace("SampleName=","SampleName=flips_data")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=control_2los1tau")
                    elif 'obs_' in ch:
                        pset=pset.replace("SampleName=","SampleName=data_obs")
                        pset=pset.replace("SelectionRegion=",
                                          "SelectionRegion=signal_2lss1tau")
                    else:
                        print 'WARNING: channel name is not valid!!!'
                        exit()

                    if '_2016g' in ch or '_2016h' in ch:
                        pset=pset.replace("HIPSafeMediumMuon=True","HIPSafeMediumMuon=False")
                    if '_2016h' in ch:
                        pset=pset.replace("Is2016H=False","Is2016H=True")
                        
                    break
                         
    vd = locals()
    selection_str = '' if args.selection=='' else '_'+args.selection
    vd['name'] = ch+'_'+args.analysis+selection_str
    vd['dataset'] = sample
    vd['cfgparams'] = pset
    vd['unit'] = perjob
    vd['outdir'] = args.outdir
    if 'data' in ch:
        vd['lumimask'] = "config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'"
        vd['runrange'] = 'config.Data.runRange = '+ run
        vd['splitting'] = 'LumiBased'
    else:
        vd['lumimask'] = ''
        vd['runrange'] = ''
        vd['splitting'] = 'EventAwareLumiBased'

    # write crab config files
    open('crab/crabConfig_'+vd['name']+'.py','wt').write(string % vd)

    # submit crab jobs
    if not args.dryrun:
        os.system('crab submit -c crab/crabConfig_'+vd['name']+'.py')
