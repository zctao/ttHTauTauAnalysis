#!/usr/bin/env python

import argparse
from ROOT import TFile, TTree

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help="Input event ntuple")

evgroup = parser.add_mutually_exclusive_group(required=True)
evgroup.add_argument('--events', type=str, nargs='+', help="event numbers")
evgroup.add_argument('--event_list', type=str,
                     help="Text file with list of event numbers")

parser.add_argument('-t','--treename', type=str, default="ttHtaus/eventTree")
parser.add_argument('-o','--output', type=str, default="evNtuple_pick.root")
args = parser.parse_args()

def getEventList(flist):
    f=open(flist,'r')
    tmp.f.read()
    return tmp.split('\n')

eventList = args.events if args.events is not None else getEventList(args.event_list)
print eventList

infile = TFile(args.input, 'read')
intree = infile.Get(args.treename)

outfile = TFile(args.output, 'recreate')
outtree = intree.CloneTree(0)

for ev in intree:
    if str(ev.nEvent) in eventList:
        outtree.Fill()

outtree.AutoSave()

print "output: ", args.output
