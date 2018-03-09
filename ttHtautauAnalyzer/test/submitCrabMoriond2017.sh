#!/bin/bash

# lumi mask
lumimask='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

samplelist='../data/SampleList_Moriond17.txt'

# 2lss1tau
# MC
submitCrabJobs.py $samplelist 2lss1tau signal --channel_list channels_moriond17_mc.txt --prefix m17_ --systematics jesup jesdown tesup tesdown -d
# Collision data
submitCrabJobs.py $samplelist 2lss1tau fakes --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d
submitCrabJobs.py $samplelist 2lss1tau flips --channel_list channels_moriond17_data.txt --prefix m17_  --lumimask $lumimask -d
# open the box
submitCrabJobs.py $samplelist 2lss1tau signal --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d

# 1l2tau
# MC
submitCrabJobs.py $samplelist 1l2tau signal --channel_list channels_moriond17_mc.txt --prefix m17_ --systematics jesup jesdown tesup tesdown -d
# Collision data
submitCrabJobs.py $samplelist 1l2tau fakes --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d
# open the box
submitCrabJobs.py $samplelist 1l2tau signal --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d

# 3l1tau
# MC
submitCrabJobs.py $samplelist 3l1tau signal --channel_list channels_moriond17_mc.txt --prefix m17_ --systematics jesup jesdown tesup tesdown -d
# Collision data
submitCrabJobs.py $samplelist 3l1tau fakes --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d
# open the box
submitCrabJobs.py $samplelist 3l1tau signal --channel_list channels_moriond17_data.txt --prefix m17_ --lumimask $lumimask -d
