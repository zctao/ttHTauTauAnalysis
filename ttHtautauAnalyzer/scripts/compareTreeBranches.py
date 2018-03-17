#!/usr/bin/env python

from ROOT import TFile, TCanvas, TLegend, gDirectory, gStyle, kRed

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="Input file")
parser.add_argument("branch1", type=str)
parser.add_argument("branch2", type=str)
parser.add_argument("--cut1", type=str, default="")
parser.add_argument("--cut2", type=str, default="")
parser.add_argument("--drawoption1", type=str, default="")
parser.add_argument("--drawoption2", type=str, default="")
parser.add_argument("--title", type=str)
parser.add_argument("-t","--treename", type=str, default="mva")
parser.add_argument("-o","--output", type=str, default="compare.pdf")

args = parser.parse_args()

fin = TFile(args.input, "read")
tree = fin.Get(args.treename)

tree.Draw(args.branch1+">>htemp1", args.cut1)
h1 = gDirectory.Get("htemp1")
tree.Draw(args.branch2+">>htemp2", args.cut2)
h2 = gDirectory.Get("htemp2")

#leg = TLegend()

canvas = TCanvas()

ymax = max(h1.GetMaximum(), h2.GetMaximum())
h1.SetMaximum(ymax*1.2)
h2.SetLineColor(kRed)

gStyle.SetOptStat(0)

#leg.AddEntry(h1, args.branch1+"("+args.cut1+")", "l")
#leg.AddEntry(h2, args.branch2+"("+args.cut2+")", "l")

h1.Draw(args.drawoption1)
option2 = "same" if args.drawoption2=="" else "same "+args.drawoption2
h2.Draw(option2)
#leg.Draw("same")

canvas.BuildLegend()

if args.title is None:
    t = args.input.split("/")[-1]
    t = t.replace(".root","")
    t = t.replace("mvaVars_","")
    t = t.replace("ntuple_","")
    h1.SetTitle(t)
else:
    h1.SetTitle(args.title)

canvas.SaveAs(args.output)
