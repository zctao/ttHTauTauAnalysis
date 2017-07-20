# submit multiple crab jobs
# modified from Jorge's script

import os, sys
import crab_channels

channels = crab_channels.channels_train
#samplelist = "../data/SampleList_Moriond17.txt"
samplelist = "../data/SampleList_Training.txt"

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
config.Data.outLFNDirBase = '/store/user/ztao/ttH_80X'
config.Site.storageSite = 'T3_US_FNALLPC'
'''

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
                    pset = f.next().strip()
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
                        
                    perjob = f.next().strip()
                    break
                
            else:
                if line.strip()==remove_prefix(ch,"data_obs") or line.strip()==remove_prefix(ch,"fakes_data") or line.strip()==remove_prefix(ch,"flips_data"):
                    sample = f.next().strip()
                    run = f.next().strip()
                    perjob = f.next().strip()
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
    vd['name'] = ch
    vd['dataset'] = sample
    vd['cfgparams'] = pset
    vd['unit'] = perjob
    if 'data' in ch:
        vd['lumimask'] = "config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'"
        vd['runrange'] = 'config.Data.runRange = '+ run
        vd['splitting'] = 'LumiBased'
    else:
        vd['lumimask'] = ''
        vd['runrange'] = ''
        vd['splitting'] = 'EventAwareLumiBased'

    open('crab/crabConfig_'+ch+'.py','wt').write(string % vd)
    os.system('crab submit -c crab/crabConfig_'+ch+'.py')
