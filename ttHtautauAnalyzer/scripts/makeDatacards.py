#!/usr/bin/env python

from ROOT import TFile, TH1D
import ttHTauTauAnalysis.ttHtautauAnalyzer.datacards as dc
from ttHTauTauAnalysis.ttHtautauAnalyzer.cross_section import CrossSection
import copy
import argparse

def getHistogramNames_mc(anaType, sample, addSyst, syst_coname='_CMS_ttHl_'):

    lowercase=True
    s = sample.lower() if lowercase else sample
    
    namelist = [s, s+'_gentau', s+'_faketau']
    if sample=='ttH':
        namelist += [s+'_htt', s+'_htt_gentau', s+'_htt_faketau',
                     s+'_hww', s+'_hww_gentau', s+'_hww_faketau',
                     s+'_hzz', s+'_hzz_gentau', s+'_hzz_faketau']

    if not addSyst:
        return ['x_'+hname for hname in namelist]

    systlist = []
    #systlist = ['JESUp','JESDown','TESUp','TESDown']

    for btsys in dc.BTagSysts:
        systlist.append('btag_'+btsys)
    for thu in dc.ThSysts:
        systlist.append('thu_shape_'+thu)
    for frjt in dc.FakeTauSysts:
        systlist.append(frjt)

    return ['x_'+hname for hname in namelist]+['x_'+hname+syst_coname+syst for hname in namelist for syst in systlist if not(syst in dc.FakeTauSysts and (anaType=='1l2tau' or 'gentau' in hname))]

def getHistogramNames_data(anaType, channel, addSyst, syst_coname='_CMS_ttHl_'):
    namelist = ['x_'+channel]

    if not (addSyst and 'fakes' in channel):
        return namelist

    # fake rate systematics
    systlist = []
    if anaType=='1l2tau':
        systlist = [frjt for frjt in dc.FakeTauSysts]
        # lepton fake rate systematics?
    elif anaType=='2lss1tau' or anaType=='3l1tau':
        systlist = [frl for frl in dc.FakeRateLepSysts]
        # closure test
        #systlist += [clos for clos in dc.ClosureTests]
        # closure test shapes and dealt with separately

    return namelist+[hname+syst_coname+syst for hname in namelist for syst in systlist]

def getShapefromSample_mc(anaType, sample, histname, treename, nbin, xmin, xmax,
                          ntuplelist, binningMap, lumi):
    # determine jet/tau energy correction
    correction = None
    for cor in ['JESUp','JESDown','TESUp','TESDown']:
        if cor in histname:
            correction = cor
            break
   
    # open ntuple file and get tree
    infile = dc.getNtupleFileName_mc(ntuplelist, anaType, sample, correction)
    f = TFile(infile,'read')
    tree = f.Get(treename)

    inverseSumGenWeight = tree.GetWeight()
    xsection = CrossSection[sample]

    # make data cards
    shape = dc.getShapeFromTree(tree, histname, nbin, xmin, xmax, binningMap)
    # scale histogram
    shape.Scale(lumi*xsection*inverseSumGenWeight)
    dc.makeBinContentsPositive(shape)

    shape.SetDirectory(0)
    
    return shape

def getShapefromSamples_data(anaType, channel, histname, treename, nbin, xmin, xmax,
                            ntuplelist, binningMap, lumi):
    # open ntuple files and get trees
    inputfiles = [TFile(dc.getNtupleFileName_data(ntuplelist, anaType, channel, sample),'read') for sample in dc.SamplesInChannel[channel]]
    trees = [fin.Get(treename) for fin in inputfiles]

    shape = dc.getShapeFromMergingTrees(trees, histname, nbin,xmin,xmax, binningMap)

    # make bin contents positive
    dc.makeBinContentsPositive(shape)

    shape.SetDirectory(0)
    return shape
    


            

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('anaType', choices=['1l2tau','2lss1tau','3l1tau'],
                        help="Analysis type")
    parser.add_argument('channels', nargs='+',choices=['ttH','TTW','TTZ','EWK','Rares','fakes_data','flips_data','data_obs'],
                        help="List of channels to make data cards")
    parser.add_argument('nbin', type=int, help="number of bins")
    parser.add_argument('--xmin', type=float)
    parser.add_argument('--xmax', type=float)
    parser.add_argument('-n','--ntuplelist', type=str, default='mvaNtuplelist.txt')
    parser.add_argument('-b','--binmap', type=str, default='binning.root',
                        help="Binning Map from 2D to 1D. Default histogram name: 'hTargetBinning'")
    parser.add_argument('-o','--outname', type=str, default='./datacards.root',
                        help="Output file name")
    parser.add_argument('-l', '--luminosity', type=float, default=100.,
                        help="Integrated luminosity")
    parser.add_argument('--treename', type=str, default="mva",
                        help="Name of tree")
    parser.add_argument('--sys_coname', type=str, default='_CMS_ttHl_',
                        help="common string in histogram names for systematics")
    parser.add_argument('-s','--systematics', action='store_true',
                        help="Include systematics")
    parser.add_argument('-v', '--verbose', action='count', #action='store_true',
                        help='verbosity')
    
    args = parser.parse_args()

    # get binning Map 2D histogram
    fMap = TFile(args.binmap)
    hBinningMap = fMap.Get("hTargetBinning")

    # determine xmin and xmax if not specified
    xMin = 0. if args.xmin is None else args.xmin
    xMax = float(args.nbin) if args.xmax is None else args.xmax
    
    datacards = []
    
    for channel in args.channels:
        
        print "======================================="
        print "processing channel:", channel
        print "---------------------------------------"

        
        if 'data' in channel: # Collision data

            hnames = getHistogramNames_data(args.anaType, channel, args.systematics)
            for histname in hnames:

                h = getShapefromSamples_data(args.anaType, channel, histname,
                                             args.treename, args.nbin, xMin, xMax,
                                             args.ntuplelist, hBinningMap,
                                             args.luminosity)
                datacards.append(h)

                if args.verbose>=1 and histname=='x_'+channel:
                    yields = h.Integral()
                    print channel, "\t", '%.5f'%yields
                
                # Closure test
                if not (histname=="x_fakes_data" and args.systematics):
                    continue

                if not (args.anaType=="2lss1tau" or args.anaType=="3l1tau"):
                    continue
                
                # FIXME get auxiliary files
                fname_closure = "../data/Closure_FR_syst/Closure_FR_lepton_syst_2lSS1tau_nofaketau_MVA_2lSS_36.8fb.root" if args.anaType=='2lss1tau' else "../data/Closure_FR_syst/Closure_FR_lepton_syst_3l1tau_nofaketau_MVA_3l_36.8fb.root"
                
                for clos in dc.ClosureTests:
                    h_clos = dc.getClosureTestShape(h, clos, fname_closure)
                    datacards.append(h_clos)
                     
        else: # Monte Carlo
            
            datacards_channel = []
            samples = dc.SamplesInChannel[channel]

            if args.verbose>=1:
                print '\tyields\tyields(gentau)\tyields(faketau)'
            
            firstsample = True
            for sample in samples:

                datacards_sample = []
                
                hnames = getHistogramNames_mc(args.anaType, sample, args.systematics)
                for histname in hnames: 
                    h = getShapefromSample_mc(args.anaType, sample, histname,
                                              args.treename, args.nbin, xMin, xMax,
                                              args.ntuplelist, hBinningMap,
                                              args.luminosity)
                    datacards_sample.append(h)

                if args.verbose>=2 and len(samples)>1:
                    dc.printYields(sample.lower(), datacards_sample)

                if firstsample:
                    datacards_channel = [h.Clone(h.GetName().replace(sample.lower(),channel)) for h in datacards_sample]
                else:
                    for hc, hs in zip(datacards_channel, datacards_sample):
                        hc.Add(hs)
                        
                firstsample = False

                if len(samples)>1:
                    datacards += datacards_sample   
                
            assert(datacards_channel!=[])
            datacards += datacards_channel

            if args.verbose>=1:
                if channel=='ttH':
                    dc.printYields('ttH_htt', datacards_channel)
                    dc.printYields('ttH_hww', datacards_channel)
                    dc.printYields('ttH_hzz', datacards_channel)

                dc.printYields(channel, datacards_channel)
                
                
    # write datacards to file
    outfile = TFile(args.outname, 'recreate')

    for datacard in datacards:
        datacard.Write()

    print "Output file:", args.outname
