#!/usr/bin/env python

from ROOT import TFile, TH1D
from ttHTauTauAnalysis.ttHtautauAnalyzer.plotting import DrawStackHistograms
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="datacards finle name")
parser.add_argument("-o", "--outname", type=str, default="hstack.pdf",
                    help="output plot name")
parser.add_argument("-c","--channels", nargs='+',
                    choices=['ttH','TTW','TTZ','EWK','Rares','fakes_data','data_obs'],
                    default=['Rares','EWK','TTW','TTZ','ttH','fakes_data','data_obs'])
parser.add_argument("-u", "--unblind", action='store_true')
parser.add_argument("-v", "--verbose", action='store_true')

args = parser.parse_args()

histograms = []

fin = TFile(args.input, "read")
for ch in args.channels:
    h = fin.Get('x_'+ch)
    histograms.append(h)

DrawStackHistograms(histograms, args.outname, not args.unblind, args.verbose)
