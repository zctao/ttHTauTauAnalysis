#!/usr/bin/env python

from ROOT import TFile
from ttHTauTauAnalysis.ttHtautauAnalyzer.plotting import DrawRatioStackHistograms
import ttHTauTauAnalysis.ttHtautauAnalyzer.datacards as dc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("variable", type=str, help="Variable to plot")
parser.add_argument("infile", type=str, help="Input file name")
parser.add_argument("-t","--title", type=str, help="plot title")
parser.add_argument("-o", "--outname", type=str, default="plot.pdf",
                    help="Output plot name")
parser.add_argument("-c","--channels", nargs='+', choices=dc.Channels)
parser.add_argument("-v", "--verbose", action='store_true')
parser.add_argument("-x", "--xtitle", type=str, default="BDT")
parser.add_argument("-y", "--ytitle", type=str, default="Events")

args = parser.parse_args()

def SetHistogramStyle(h, channel):
    h.SetLineColor(1)
    h.SetLineWidth(1)
    
    if channel=="data_obs":
        h.SetMarkerColor(1)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(2)
        h.SetLineStyle(1)
    elif channel=="ttH":
        h.SetFillColor(628)
    elif channel=="TTZ":
        h.SetFillColor(822)
    elif channel=="TTW":
        h.SetFillColor(823)
    elif channel=="TTWW":
        h.SetFillColor(824)
    elif channel=="EWK":
        h.SetFillColor(610)
    elif channel=="Rares":
        h.SetFillColor(851)
    elif channel=="fakes_data":
        h.SetFillColor(1)
        h.SetFillStyle(3005)
    elif channel=="flips_data":
        h.SetFillColor(1)
        h.SetFillStyle(3004)
    elif channel=="conversions" or channel=="Convs":
        h.SetFillColor(5)
        

histogram_obs = None
histograms_pred = []

fin = TFile.Open(args.infile)
for key in fin.GetListOfKeys():
    h = key.ReadObj()
    if not 'TH1' in h.ClassName():
        continue
    
    for ch in args.channels:
        if args.variable+'_'+ch == h.GetName():

            SetHistogramStyle(h,ch)
            
            if ch == "data_obs":
                histogram_obs = h
            else:
                h.SetTitle(ch)
                histograms_pred.append(h)

DrawRatioStackHistograms(args.outname, args.title, histograms_pred, histogram_obs,
                         xTitle=args.xtitle, yTitle_top=args.ytitle,
                         yTitle_bot="#frac{Data}{Expectation}",
                         ymin=0.1, ymax=10000, logScale=True)
