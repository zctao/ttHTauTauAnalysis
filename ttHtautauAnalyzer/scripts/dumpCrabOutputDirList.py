#!/usr/bin/env python

import os
import argparse

def getSampleFullname(sample, list):

    with open(list) as f:
        for line in f:
            if sample != line.strip():
                continue

            fullname = f.next().strip()
            return fullname.split('/')[1] # [0] is expected to always be "'"

        
def getSubDirectoryNames(eosDir):
    # return a list of sub-directory names under eos
    subdirs = os.popen('eos root://cmseos.fnal.gov ls ' + eosDir).read().strip().split('\n')
    return subdirs


def getLatestTimeStamp(timestamplist):

    outint = 0
    outstr = ''
    
    for ts in timestamplist:       
        ts_str = ts.replace('_','')
        assert(len(ts_str)==12)
        its = int(ts_str)

        if its > outint:
            outint = its
            outstr = ts

    return outstr
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Dump the list eos directories of crab job outputs given samples names and analysis type')

    parser.add_argument('type', type=str, choices=['1l2tau','2lss1tau','3l1tau'],
                        help='analysis type')
    parser.add_argument('samples', nargs='+', type=str,
                        help="sample names. NOTE: They will be hadd'ed all together into ONE root file. Normally they should be the same type of samples e.g. TTW and TTW_ext")
    parser.add_argument('-l','--list', type=str, default='SampleList.txt',
                        help="Sample list")
    parser.add_argument('-o','--output', type=str, default='list.txt',
                        help="Output name")
    parser.add_argument('-t','--timestamp', type=str, default='',
                        help="Time stamp of crab job. If none provided, assume the latest in the directory")
    parser.add_argument('-r','--rootdir', type=str,
                        default='/store/user/ztao/ttH_80X/',
                        help="root directory in eos")
    parser.add_argument('-s','--streaming', action='store_true',
                        help="Stream out list of directories instead of write to disk")
    args = parser.parse_args()

    outfile = open(args.output,'w')
    
    for sample in args.samples:
        
        fullsamplename = getSampleFullname(sample, args.list)
        directory = args.rootdir+fullsamplename+'/crab_'+sample+'_'+args.type+'/'
        
        # time stamp
        if args.timestamp!='':
            directory += args.timestamp+'/'
        else:
            # loop over the sub directories, which should be the time stampes at this level, and pick the latest one
            TimeStampList = getSubDirectoryNames(directory)
            directory += getLatestTimeStamp(TimeStampList)+'/'
                
        # 0000, 0001,...
        sublist = getSubDirectoryNames(directory)      
        
        for number in sublist:
            outstr = directory + number + '/'
            
            if args.streaming:
                print outstr
            else:
                outfile.write(outstr+'\n')
                
    outfile.close()

    if args.streaming:
        os.remove(args.output) # since it's empty
