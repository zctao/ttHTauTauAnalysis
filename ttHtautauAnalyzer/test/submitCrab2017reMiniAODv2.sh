#!/bin/bash

# lumimask
lumimask='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
# sample list
datasets='../dataFiles/DatasetList_2017reMiniAODv2.csv'

# MC
#submitCrabJobsv2.py $datasets --prefix 2018jun_ --samples ttHJetToNonbb TTZ TTGJets TTW TTW_psw WZ TTToDiLep TTToDiLep_psw TTToSemiLep TTToSemiLep_psw TTToHad TTToHad_psw #ZZ ZZ_ext WW WWW WWZ WZZ ZZZ TTTT tZq WWds -d

# Data
#submitCrabJobsv2.py $datasets --prefix 2018jun_ --sample_list samples_2017remaodv2_data.txt --lumimask $lumimask -d
#submitCrabJobsv2.py $datasets --prefix 2018jun_ --samples data_e_2017b data_e_2017c data_e_2017d data_e_2017e data_e_2017f --lumimask $lumimask -d
