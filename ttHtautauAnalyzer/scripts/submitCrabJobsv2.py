#!/usr/bin/env python

import os
import argparse

from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import channel_suffix, getDatasetDict, getSamples

parser = argparse.ArgumentParser()

parser.add_argument('DatasetList', type=str, help="List of Datasets")

chgroup = parser.add_mutually_exclusive_group(required=True)
chgroup.add_argument('--samples',nargs='+', type=str,
                    help="samples to run crab jobs on")
chgroup.add_argument('--sample_list', type=str,
                     help="File that stores a list of samples to run crab jobs on")

parser.add_argument('--analysis_label', choices=['incl','gen'], default='incl',
                    help="incl for ntuple production. gen for genParticle analysis")
parser.add_argument('--prefix', type=str, default='',
                    help="Prefix of crab job name.")
parser.add_argument('-c','--config', type=str, default='analyzer2017_cfg.py',
                    help="cmsRun config file")
parser.add_argument('--dbs', choices=['global','phys01','phys02','phys03'],
                    default='global',
                    help="Alias of the URL of the DBS reader instance where the input dataset is published ")
#parser.add_argument('--systematics', nargs='+',
#                    choices=['nominal','jesup','jesdown','tesup','tesdown'],
#                    default=['nominal'], help="Energy correction")
parser.add_argument('-o','--outdir', type=str,
                    default='/store/user/ztao/ttHtaus_94X', help="Output directory")
parser.add_argument('-d','--dryrun', action='store_true',
                    help="Run the script without actually submitting crab jobs")
parser.add_argument('--lumimask', type=str,
                    default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt', help='Lumi mask to run on data sample')
parser.add_argument('-a','--automatic', action='store_true',
                    help="Automatic data splitting")
parser.add_argument('-b','--batch_lpc', action='store_true',
                    help="Submit batch job on CMS LPC")

args = parser.parse_args()

def getParametersetString(sample, isdata):
    pset = "['SampleName=','isData=False','doCutFlow=True']"
    
    # sample name
    sname = sample
    for suffix in channel_suffix:
        sname=sname.replace(suffix,'')
    pset=pset.replace("SampleName=","SampleName="+sname)
    
    # is collision data
    if isdata:
        pset=pset.replace("isData=False","isData=True")

    return pset

dataset_dict = getDatasetDict(args.DatasetList)
samples = args.samples if args.samples is not None else getSamples(args.sample_list)


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
config.JobType.inputFiles=['qg']
config.Data.inputDataset = %(dataset)s
config.Data.inputDBS = '%(DBS)s'
config.Data.splitting = '%(splitting)s'
config.Data.unitsPerJob = %(unit)s
%(lumimask)s
%(runrange)s
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = '%(outdir)s'
config.Site.storageSite = 'T3_US_FNALLPC'
'''

if args.automatic:
    template = template.replace("config.Data.unitsPerJob","#config.Data.unitsPerJob")

if args.batch_lpc:
    template = template.replace("config.Data.ignoreLocality = False",
                                "config.Data.ignoreLocality = True\n"
                                +"config.Site.whitelist = ['T3_US_FNALLPC']\n"
                                +"config.Site.ignoreGlobalBlacklist = True\n")

for sample in samples:
    print sample
    isData = 'data' in sample
        
    jobname = args.prefix + '_' + sample + '_'+args.analysis_label

    parameterset = getParametersetString(sample, isData)
            
    vd = locals()
    vd['name'] = jobname
    vd['dataset'] = "'"+dataset_dict[sample]['dataset']+"'"
    vd['configfile'] = args.config
    vd['cfgparams'] = parameterset
    vd['DBS'] = args.dbs
    vd['unit'] = dataset_dict[sample]['unitPerJob']
    vd['outdir'] = args.outdir
    if isData:
        vd['lumimask'] = "config.Data.lumiMask = "+"'"+args.lumimask+"'"
        vd['runrange'] = 'config.Data.runRange = '+ "'"+dataset_dict[sample]['runRange']+"'"
        vd['splitting'] = 'LumiBased'
    else:
        vd['lumimask'] = ''
        vd['runrange'] = ''
        vd['splitting'] = 'EventAwareLumiBased'
        
    if args.automatic:
        vd['splitting'] = 'Automatic'

    # write crab config files
    open('crab/crabConfig_'+vd['name']+'.py','wt').write(template % vd)

    # submit crab jobs
    if not args.dryrun:
        os.system('crab submit -c crab/crabConfig_'+vd['name']+'.py')
