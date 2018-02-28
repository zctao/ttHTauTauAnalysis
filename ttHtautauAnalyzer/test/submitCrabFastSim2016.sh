#!/bin/bash

samplelist='../data/SampleList_FastSim2016.txt'

# 2lss1tau
submitCrabJobs.py $samplelist 2lss1tau loose --channel_list channels_fastsim16.txt --prefix fs16_ --dbs phys03 -a -b -d

# 1l2tau
submitCrabJobs.py $samplelist 1l2tau loose --channel_list channels_fastsim16.txt --prefix fs16_ --dbs phys03 -a -b -d

# 3l1tau
submitCrabJobs.py $samplelist 3l1tau loose --channel_list channels_fastsim16.txt --prefix fs16_ --dbs phys03 -a -b -d
