#!/usr/bin/env python
from ROOT import TFile, TH1D, gROOT
gROOT.SetBatch(True)
import ttHTauTauAnalysis.ttHtautauAnalyzer.datacards as dc
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
parser.add_argument('--sysname', type=str, default='_CMS_ttHl_',
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

    xbins = dc.getBinEdges(args.analysis)
    ngroup = len(xbins)-1
    hname = hist.GetName()

    # rename if necessary
    hname = updateHistName(hname)
    
    #hist.Sumw2()
    hist_rebinned = hist.Rebin(ngroup, hname, xbins)
    
    if args.positvebins:
        dc.makeBinContentsPositive(hist_rebinned, verbosity=args.verbose)

    # crop all uncertainties to 100% to avoid negative variations
    for bin in xrange(1, hist_rebinned.GetNbinsX()+1):
        if hist_rebinned.GetBinError(bin) > hist_rebinned.GetBinContent(bin):
            hist_rebinned.SetBinError(bin, hist_rebinned.GetBinContent(bin))
    
    if args.scalebywidth:
        hist_rebinned.Scale(1.,'width')
        # set y axis label?

    histogramsRebinned.append(hist_rebinned)

# write to file
outfile = TFile(args.output, 'recreate')
for h in histogramsRebinned:
    h.Write()

print "Output file:", args.output
