#!/usr/bin/env python

import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('label', type=str, help="Crab job label for the input event ntuples. E.g. 2019jan")
parser.add_argument('samples', nargs='+', help="sample names")
parser.add_argument('--datasetlist', type=str,
                    default="DatasetList_2017reMiniAODv2.csv")
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
parser.add_argument('-c','--corrections', nargs='+',
                    choices=['NA','JESUp','JESDown','JERUp','JERDown',
                             'TESUp','TESDown'],
                    default=['NA'], help="Jet/Tau energy correction")
parser.add_argument('-m','--genmatching', action='store_true',
                    help="Require gen matching in event selections")
parser.add_argument('-l', '--log', type=str, help="Log name")
parser.add_argument('-t','--topeosdir', type=str,
                    default="/store/user/ztao/ttHtaus_94X/")
parser.add_argument('-o','--outdir', type=str, help="Output directory")

parser.add_argument('--condor_storage', type=str,
                    default='/store/user/ztao/Condor/mvaNtuples/')
parser.add_argument('--tar_dir', type=str, default='/store/user/ztao/Condor/',
                    help="Directory in which the cmssw tar file is stored")
parser.add_argument('--cmssw', type=str, default='CMSSW_9_4_7',
                    help="CMSSW version")
parser.add_argument('--arch', type=str, default='slc6_amd64_gcc630',
                    help='scram arch')
parser.add_argument('--fname_jdl', type=str, default='condor_mvantuple.jdl')
parser.add_argument('--fname_sh', type=str, default='makeMVANtuples_condor.sh')

args = parser.parse_args()

run_command = "produceMVANtuplesv2.py"
run_command += ' '+args.label
for sample in args.samples:
    run_command += ' '+sample
run_command += ' --datasetlist '+args.cmssw+'/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/'+args.datasetlist

run_command += ' --analysis'
for ana in args.analysis:
    run_command += ' '+ana

run_command += ' --jobtypes'
for job in args.jobtypes:
    run_command += ' '+job

if args.selection is not None:
    run_command += ' --selection '+args.selection

run_command += ' -c'
for correction in args.corrections:
    run_command += ' '+correction

run_command += ' --topeosdir '+args.topeosdir
    
if args.genmatching:
    run_command += ' -m'

if args.log is not None:
    run_command += ' --log '+args.log

outntuple_dir = args.label if args.outdir is None else args.outdir
run_command += ' -o '+outntuple_dir

run_command += ' --transfer_inputs'

print run_command

# Condor Job Description template
template_jdl = '''universe = vanilla
Executable = %(exe)s
Error = log/job.error.$(Cluster)-$(Process)
Output = log/job.output.$(Cluster)-$(Process)
Log = log/job.log.$(Cluster)-$(Process)
stream_output = false
stream_error = false
notification = never
Transfer_Input_Files = %(inputs)s
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
'''

# Executable template
template_sh = '''#!/bin/bash

workdir=$(pwd)
echo Work directory $workdir
cmssw=%(CMSSW)s
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=%(scram_arch)s

xrdcp -s root://cmseos.fnal.gov/%(tardir)s/$cmssw.tgz .
tar -xf $cmssw.tgz
rm $cmssw.tgz
cd $cmssw/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`
scramexit=$?
if [[ $scramexit -ne 0 ]]; then
    echo scram exited with code $scramexit
    exit $scramexit
fi

echo "CMSSW: "$CMSSW_BASE

cd $workdir
%(run_command)s

echo "*******************************************"
# copy output to eos
# all output root files should be under $workdir/%(ntupledir)s/
#eos root://cmseos.fnal.gov mkdir -p %(fulloutdir)s
outeosdir=root://cmseos.fnal.gov/%(fulloutdir)s
cd $workdir/%(ntupledir)s/
for file in *.root
do
        echo "xrdcp -f ${file} ${outeosdir}/${file}"
        xrdcp -f ${file} ${outeosdir}/${file}
        xrdexit=$?
        if [[ $xrdexit -ne 0 ]]; then
                rm *.root
                echo xrdcp exited with code $xrdexit
                exit $xrdexit
        fi
        rm ${file}      
done
cd $workdir
rm -rf $cmssw
rm -rf %(ntupledir)s/
'''



# jdl file
vd = locals()
vd['exe'] = args.fname_sh
vd['inputs'] = args.fname_sh

# bash script
vd = locals()
vd['CMSSW'] =  args.cmssw
vd['scram_arch'] = args.arch
#fname_tar = args.tar_dir.rstrip('/')+'/'+args.cmssw+'.tgz'
#vd['tarfile'] = fname_tar
vd['tardir'] = args.tar_dir.rstrip('/')
vd['run_command'] = run_command
vd['fulloutdir'] = args.condor_storage.rstrip('/')+'/'+outntuple_dir
vd['ntupledir'] = outntuple_dir

# write jdl and bash script
open(args.fname_jdl,'wt').write(template_jdl % vd)
print "Generated jdl file: "+args.fname_jdl
open(args.fname_sh,'wt').write(template_sh % vd)
print "Generated bash script: "+args.fname_sh
print "To submit the job: condor_submit "+args.fname_jdl
print "Make sure a tarball of the current CMSSW directory ("+args.cmssw+".tgz) is in "+args.tar_dir.rstrip('/')+"/ before submitting the condor job."
print "Also make sure there is a log/ directory where the condor jobs are submitted."
