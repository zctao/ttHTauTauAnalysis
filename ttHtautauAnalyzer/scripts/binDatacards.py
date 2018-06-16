#!/usr/bin/env python
from ROOT import TFile, TH1D, gROOT
gROOT.SetBatch(True)
import ttHTauTauAnalysis.ttHtautauAnalyzer.datacards as dc
import array
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('analysis', choices=['1l2tau','2lss1tau','3l1tau','2l2tau'],
                    help="Analysis type")
parser.add_argument('inputfile', type=str, help="Input file")
parser.add_argument('-o','--output', type=str, default="datacards_binned.root",
                    help="Output name")
parser.add_argument('-c','--channels', nargs='+', choices=dc.Channels,
                    help="Channels to bin datacards")
parser.add_argument('-p','--positvebins', action='store_true',
                    help="Make bin contents positive")
parser.add_argument('-w','--scalebywidth', action='store_true',
                    help="Scale histograms by bin width")
parser.add_argument('--nosystematics', action='store_true',
                    help="Only include nominal ones")
parser.add_argument('--sysname', type=str, default='_CMS_ttHl17_',
                    help="label for systematics")
parser.add_argument('-v','--verbose', action='store_true')

args = parser.parse_args()

histogramsToRebin=[]        

# open root file and get histograms
infile = TFile.Open(args.inputfile)
for key in infile.GetListOfKeys():
    h = key.ReadObj()
    
    if not 'TH1' in h.ClassName():
        continue

    if args.nosystematics and args.sysname in h.GetName():
        continue
    
    for channel in args.channels:
        if 'x_'+channel in h.GetName():

            skip = False
            if args.sysname in h.GetName(): # for shape systematics
                if '_thu_shape_' in h.GetName() and not (channel in ['ttH','TTW','TTZ']):
                    skip = True
            if not skip:
                histogramsToRebin.append(h)
                
            break

def getBinEdges(anatype):
    binEdges=[]
    if anatype=='1l2tau':
        binEdges = dc.getUniformBinEdges(7, 0.0, 1.0)
    elif anatype=='2lss1tau':
        binEdges = [0.0, 0.14, 0.18, 0.22, 0.28, 0.32, 0.35, 0.38, 0.43, 0.47, 0.53, 1.0]         
    elif anatype=='3l1tau':
        binEdges = [0.0, 0.28, 0.35, 0.40, 0.47, 0.53, 1.0]
    elif anatype=='2l2tau':
        binEdges = [0.0, 0.35, 0.41, 0.47, 1.0]
    else:
        print "Unknow analysis type!"

    return array.array('d', binEdges)

def updateHistName(hname):
    newname = hname
    
    if '_TES' in hname:
        newname = hname.replace('_TES','_tauES')
    if 'x_ttH' in hname and '_thu_shape_' in hname:
        newname = hname.replace('_thu_shape_', '_thu_shape_ttH_')
        newname = newname.replace(args.sysname,'_CMS_ttHl_')
    if 'x_TTW' in hname and '_thu_shape_' in hname:
        newname = hname.replace('_thu_shape_', '_thu_shape_ttW_')
        newname = newname.replace(args.sysname,'_CMS_ttHl_')
    if 'x_TTZ' in hname and '_thu_shape_' in hname:
        newname = hname.replace('_thu_shape_', '_thu_shape_ttZ_')
        newname = newname.replace(args.sysname,'_CMS_ttHl_')

    return newname

histogramsRebinned=[]

for hist in histogramsToRebin:
    print hist.GetName()

    xbins = getBinEdges(args.analysis)
    ngroup = len(xbins)-1
    hname = hist.GetName()

    # rename if necessary
    hname = updateHistName(hname)
    
    #hist.Sumw2()
    hist_rebinned = hist.Rebin(ngroup, hname, xbins)
    
    if args.positvebins:
        dc.makeBinContentsPositive(hist_rebinned, verbosity=args.verbose)    

    if args.scalebywidth:
        hist_rebinned.Scale(1.,'width')
        # set y axis label?
        
    histogramsRebinned.append(hist_rebinned)

# write to file
outfile = TFile(args.output, 'recreate')
for h in histogramsRebinned:
    h.Write()

print "Output file:", args.output
