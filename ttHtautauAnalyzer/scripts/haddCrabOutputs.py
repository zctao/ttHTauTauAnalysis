#!/usr/bin/env python

import os, subprocess
import argparse
from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname, getSubDirectoryNames

def getLatestTimeStamp(timestamplist, maxtimestamp='999999_999999'):

    outint = 0
    outstr = ''

    maxts_str = maxtimestamp.replace('_','')
    assert(len(maxts_str)==12)
    maxts = int(maxts_str)
    
    for ts in timestamplist:       
        ts_str = ts.replace('_','')
        assert(len(ts_str)==12)
        its = int(ts_str)

        if its > outint and its < maxts:
            outint = its
            outstr = ts

    return outstr

def getCrabOutputDirectories(sample, datasetList, crabLabel,
                             jobType='incl', timestamp='',
                             eosDir='/store/user/ztao/ttHtaus_94X'):
    outdir = []
    
    # EOS top directory
    directory = eosDir.rstrip('/')+'/'
    
    # Full name of the sample as in CMSDAS
    directory += getSampleFullname(sample, datasetList)+'/'
    
    # Crab job directory
    directory += 'crab_'+crabLabel+'_'+sample+'_'+jobType+'/'

    # Time stamp
    # Loop over the sub directories, which should be the time stampes of the jobs
    TimeStampList = getSubDirectoryNames(directory)
    assert(len(TimeStampList)>=1)
    if TimeStampList == ['']:
        raise ValueError(directory+" is empty or does not exist!!")
    #
    if timestamp == '': # Nothing specified. Get the latest one.
        directory += getLatestTimeStamp(TimeStampList)+'/'
    elif timestamp in TimeStampList: # just use it
        directory += timestamp+'/'
    else: # cannot find the specified time stamp
        raise ValueError("Cannot find "+timestamp+" in "+directory)

    # 0000, 0001, ...
    sublist = getSubDirectoryNames(directory)

    for number in sublist:
        outstr = directory + number + '/'
        outdir.append(outstr)
            
    return outdir
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('datasetList', type=str,
                        help="The CSV file containing a list of datasets")
    parser.add_argument('crabJobLabel', type=str,
                        help="Label of the crab job. Example: 2019jan. The crab project directory is expected to have the name 'crab_<crabJobLabel>_<sample>_<type>/'.")
    parser.add_argument('samples', nargs='+', type=str,
                        help="sample names. NOTE: They will be hadd'ed all together into ONE root file. Normally they should be the same type of samples e.g. TTW and TTW_ext")
    parser.add_argument('-o','--outputlabel', type=str, default=None)
    parser.add_argument('-e','--eosdir', type=str,
                        default="/store/user/ztao/ttHtaus_94X",
                        help="EOS directory")
    parser.add_argument('-a','--analysis_label', choices=['incl','gen'],
                        default='incl',
                        help="incl for ntuple production. gen for genParticle analysis")
    parser.add_argument('-t','--timestamp', type=str, default='',
                        help="Time stamp of crab job. If none provided, assume the latest in the directory")
    parser.add_argument('-w','--workspace', type=str, default='/uscmst1b_scratch/lpc1/3DayLifetime/ztao/',
                        help="hadd working space")

    args = parser.parse_args()

    # Output file name
    olabel = 'data' if 'data' in args.samples[0] else args.samples[0]
    if args.outputlabel is not None:
        olabel = args.outputlabel
    outputname = 'ntuple_'+olabel+'_'+args.analysis_label

    # Loop over samples and get a list of directories where root files are stored
    outfiledir_list = []
    for sample in args.samples:
        outfiledir_list += getCrabOutputDirectories(sample, args.datasetList,
                                                    args.crabJobLabel,
                                                    args.analysis_label,
                                                    args.timestamp, args.eosdir)
    # Write the list to disk
    tmpworkdir = args.workspace.rstrip('/')+'/'+args.crabJobLabel
    subprocess.call(['mkdir', '-p',tmpworkdir])
    fname_outlist = tmpworkdir+'/'+outputname+'_hadd_directory.txt'
    f_outlist = open(fname_outlist,'w')
    for outdir in outfiledir_list:
        f_outlist.write(outdir+'\n')

    f_outlist.close()
        
    print "Run haddEOSRoot.sh to add the root files"
    fname_outroot = tmpworkdir+'/'+outputname+'.root'
    
    subprocess.call(["haddEOSRoot.sh",fname_outroot,fname_outlist])

    # Check if the root file is successfully created
    success = os.path.isfile(fname_outroot)

    if not success:
        print "No output root file is created"
        os.remove(fname_outlist) # since hadd is not successful
    else:
        print "Copy files to the EOS space"
        # Transfer the files from the hadd workspace to eos
        # create the output director
        outputdir = args.eosdir.rstrip('/')+'/'+getSampleFullname(args.samples[0],args.datasetList)+'/'+args.crabJobLabel+'/'
        os.system('eos root://cmseos.fnal.gov mkdir '+outputdir)
        # copy
        subprocess.call(['xrdcp','-f',fname_outroot,'root://cmseos.fnal.gov/'+outputdir])
        subprocess.call(['xrdcp','-f',fname_outlist,'root://cmseos.fnal.gov/'+outputdir])

        print outputdir+outputname+'.root'
