#!/bin/bash

# lumimask
lumimask='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
# sample list
datasets='../data/DatasetList_2017reMiniAODv2.csv'

# MC
submitCrabJobsv2.py $datasets --prefix 2018apr_ --samples ttHJetToNonbb TTZ TTW TTW_psw WZ TTToDiLep TTToDiLep_psw TTToSemiLep TTToSemiLep_psw TTToHad TTToHad_psw --systematics nominal jesup jesdown tesup tesdown -d

# Data
submitCrabJobsv2.py $datasets --prefix 2018apr_ --sample_list samples_2017remaodv2_data.txt --lumimask $lumimask -d
