#!/usr/bin/env python
from ROOT import TFile, TH1D, TTree, gROOT
gROOT.SetBatch(True)
from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getDatasetDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('analysis', choices=['2lss1tau','3l1tau'])
parser.add_argument('samples', nargs='+', choices=['TTToDiLep','TTToSemiLep'])
parser.add_argument('-l','--luminosity', type=float, default=41.53,
                    help="Integrated luminosity in fb^-1")
parser.add_argument('--variable', type=str, default="mvaOutput",
                    help="Variable to make the histograms")
parser.add_argument('--xmin', type=float, default=0.,
                    help="lower limit of histograms")
parser.add_argument('--xmax', type=float, default=1.,
                    help="upper limit of histograms")
parser.add_argument('--nbins', type=int, default=100, help="Number of bins")
parser.add_argument('-o','--output', type=str, default="Closure_FR_lepton_syst.root")
parser.add_argument('-d','--datasetlist', type=str, default="DatasetList.csv")
parser.add_argument('--treename', type=str, default="mva", help="Input tree name")
parser.add_argument('--prefix', type=str, default="mvaNtuple_",
                    help="Common prefix of mva ntuple names")
parser.add_argument('--version', type=str, default="jun2018v2",
                    help="mvaNtuple version")
parser.add_argument('--eostopdir', type=str,
                    default="/store/user/ztao/Condor/mvaNtuples/",
                    help="EOS directory where mva ntuples are stored")

args = parser.parse_args()

def getShapeFromTree(tree, event_weight, histname):
    h = TH1D(histname, "", args.nbins, args.xmin, args.xmax)
    h.Sumw2()

    for ev in tree:
        value = ev.GetLeaf(args.variable).GetValue()
        weight = ev.GetLeaf(event_weight).GetValue()

        # scale factors?
        scaleFactors = ['pu_weight','mc_weight','btag_sf','lepid_sf','tauid_sf',
                        'hlt_sf']
        for sf in scaleFactors:
            assert(ev.GetLeaf(sf).GetValue() > -9998.)
            weight *= ev.GetLeaf(sf).GetValue()

        h.Fill(value, weight)
        
    return h

        
# dataset list dictionary
datalist_dict = getDatasetDict(args.datasetlist)

eosDirectoryString = 'root://cmsxrootd.fnal.gov/'+args.eostopdir+args.version+'/'

firstsample=True
hist_tt = []
hist_all = []

for sample in args.samples:
    
    hist_sample = []

    # get input mva ntuple
    ntuplefilename = eosDirectoryString+args.prefix+sample+'_psw_application_fake_'+args.analysis + '.root'
    print "Open mva ntuple :", ntuplefilename
    infile = TFile.Open(ntuplefilename)
    tree = infile.Get(args.treename)
    
     # make histogram from tree
    h_fr_tt = getShapeFromTree(tree, 'event_weight_FR_TT', 'x_'+sample+'_FR_TT_MC')
    h_fr_qcd_el = getShapeFromTree(tree, 'event_weight_FR_el_QCD_mu_TT',
                                   'x_'+sample+'_FR_QCD_MC_el')
    h_fr_qcd_mu = getShapeFromTree(tree, 'event_weight_FR_el_TT_mu_QCD',
                                   'x_'+sample+'_FR_QCD_MC_mu')
    
    # scale
    SumGenWeight = tree.GetUserInfo().FindObject('SumGenWeight').GetVal()  
    xsection = tree.GetUserInfo().FindObject('xsection').GetVal() if datalist_dict is None else float(datalist_dict[sample]['xsection'])
    if xsection < 0.:
        print "WARNING: cross section for sample", sample, " is not valid!"   
    lumi = args.luminosity*1000  # convert from 1/fb to 1/pb

    h_fr_tt.SetDirectory(0)
    h_fr_qcd_el.SetDirectory(0)
    h_fr_qcd_mu.SetDirectory(0)
    
    h_fr_tt.Scale(lumi * xsection / SumGenWeight)
    h_fr_qcd_el.Scale(lumi * xsection / SumGenWeight)
    h_fr_qcd_mu.Scale(lumi * xsection / SumGenWeight)
    hist_sample.append(h_fr_tt)
    hist_sample.append(h_fr_qcd_el)
    hist_sample.append(h_fr_qcd_mu)
    
    h_fr_tt_minus_fr_qcd_el=h_fr_tt.Clone('x_'+sample+'_FR_TT_MC_minus_FR_QCD_MC_el')
    h_fr_tt_minus_fr_qcd_el.Add(h_fr_qcd_el, -1.)
    h_fr_tt_minus_fr_qcd_el.SetDirectory(0)
    hist_sample.append(h_fr_tt_minus_fr_qcd_el)
    
    h_fr_tt_minus_fr_qcd_mu=h_fr_tt.Clone('x_'+sample+'_FR_TT_MC_minus_FR_QCD_MC_mu')
    h_fr_tt_minus_fr_qcd_mu.Add(h_fr_qcd_mu, -1.)
    h_fr_tt_minus_fr_qcd_mu.SetDirectory(0)
    hist_sample.append(h_fr_tt_minus_fr_qcd_mu)
    
    if firstsample:
        hist_tt = [h.Clone(h.GetName().replace(sample,'TT')) for h in hist_sample]
        for h in hist_tt:
            h.SetDirectory(0)
    else:
        for ht, hs in zip(hist_tt, hist_sample):
            ht.Add(hs)
            
    hist_all += hist_sample
            
    firstsample=False

outfile = TFile(args.output, 'recreate')

for hist in hist_tt:
    hist.Write()

if len(args.samples)>1:
    for hist in hist_all:
        hist.Write()

print "Output file :", args.output
