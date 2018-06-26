from ROOT import THStack,TLegend,TH1,TH1D,TRatioPlot,TTree,TCanvas,TFile,TPaveText
from ROOT import gPad, gDirectory, Double
from ROOT import kRed, kBlack, kGreen, kViolet, kAzure
import array

def getLegendLabel(channel):
    label = channel    
    if channel=="ttH":
        label = "t#bar{t}H"
    elif channel=="TTZ":
        label = "t#bar{t}Z"
    elif channel=="TTW":
        label = "t#bar{t}W"
    elif channel=="TTWW":
        label = "t#bar{t}WW"
    elif channel=="EWK":
        label = "Electroweak"
    elif channel=="fakes_data":
        label = "Fakes"
    elif channel=="flips_data":
        label = "Flips"

    return label

def setUncertaintyStyle(hist, style=1001):
    alpha = 0.40
    color=12
    hist.SetLineColorAlpha(color, alpha)
    hist.SetLineWidth(0)
    hist.SetFillColorAlpha(color, alpha)
    hist.SetFillStyle(style)

def getLabels_CMS_preliminary(x0, y0, x0_luminosity, lumi="41.53"):
    label_cms = TPaveText(x0, y0 + 0.0025, x0 + 0.0950, y0 + 0.0600, "NDC")
    label_cms.AddText("CMS")
    label_cms.SetTextFont(61)
    label_cms.SetTextAlign(13)
    label_cms.SetTextSize(0.0575)
    label_cms.SetTextColor(1)
    label_cms.SetFillStyle(0)
    label_cms.SetBorderSize(0)
    #label_cms.Draw()

    label_preliminary = TPaveText(x0 + 0.120, y0 - 0.0010, x0 + 0.2950, y0 + 0.0500, "NDC");
    label_preliminary.AddText("Preliminary");
    label_preliminary.SetTextFont(52);
    label_preliminary.SetTextAlign(13);
    label_preliminary.SetTextSize(0.050);
    label_preliminary.SetTextColor(1);
    label_preliminary.SetFillStyle(0);
    label_preliminary.SetBorderSize(0);
    #label_preliminary.Draw();

    label_luminosity = TPaveText(x0_luminosity, y0 + 0.0050, x0_luminosity + 0.1900, y0 + 0.0550, "NDC");
    label_luminosity.AddText(lumi+" fb^{-1} (13 TeV)");
    label_luminosity.SetTextAlign(13);
    label_luminosity.SetTextSize(0.050);
    label_luminosity.SetTextColor(1);
    label_luminosity.SetFillStyle(0);
    label_luminosity.SetBorderSize(0);
    #label_luminosity.Draw();

    return label_cms, label_preliminary, label_luminosity


def DrawRatioStackHistograms(outname, legendTitle,
                             histograms_predicted, histogram_observe, 
                             ymin = 0.1, ymax = 10000,
                             xTitle = None, yTitle_top = None, yTitle_bot = None,
                             sortByYields=True, verbose=False, unblinded=False,
                             logScale=True, addCMSLabel=True):
    # histogram_observe: a histogram for observed data
    # histograms_predicted: a list of histograms to be plotted in THStack
    # options_observe: dictionary for drawing options
    # options_predicted: a list of dictionaries for drawing options
    
    hs = THStack("h_stack","")

    legend = TLegend(0.59,0.70,0.99,0.90)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetTextSize(0.02)
    legend.SetNColumns(2)
    legend.SetHeader(legendTitle)
    
    if histogram_observe is not None:
        legend.AddEntry(histogram_observe, "Observed", "p")

    if sortByYields:
        histograms_predicted.sort(key=lambda h: h.Integral())
    
    for hist in histograms_predicted:
        hs.Add(hist)
        #label = getLegendLabel(hist.GetName().split('_')[-1])
        label = getLegendLabel(hist.GetTitle())
        legend.AddEntry(hist, label, "f")

    h_pred = hs.GetStack().Last().Clone()
    #h_pred.SetLineColorAlpha(0,1)

    if verbose:
        pred_err = Double(0.)
        nbins_pred = h_pred.GetNbinsX()
        pred_yields = h_pred.IntegralAndError(1, nbins_pred, pred_err)
        print "Predicted yields :", pred_yields, "+/-", pred_err
        if unblinded:
            print "Observed yields :", histogram_observe.Integral()
           
    canvas = TCanvas("canvas", "canvas", 950, 1100)
    
    if logScale:
        gPad.SetLogy()

    h_pred.SetLineColorAlpha(1,0)
    h_pred.SetMarkerColorAlpha(1,0)
    h_pred.SetFillColorAlpha(1,0)
        
    if histogram_observe is None:
        histogram_observe=h_pred.Clone("dummy")
        histogram_observe.SetMarkerStyle(20)
        
    histogram_observe.SetStats(0)
    histogram_observe.GetXaxis().SetTitle(xTitle)
    histogram_observe.SetMaximum(ymax)
    histogram_observe.SetMinimum(ymin)
        
    ratioplot = TRatioPlot(histogram_observe, h_pred, "divsym")
    ratioplot.SetSeparationMargin(0.01)
    ratioplot.SetLeftMargin(0.12)
    ratioplot.SetRightMargin(0.03)
    ratioplot.SetUpTopMargin(0.12)
    #ratioplot.SetLowBottomMargin(0.30)
    ratioplot.SetH1DrawOpt("e1p")
    ratioplot.Draw()
    
    ratioplot.GetUpperRefYaxis().SetTitle(yTitle_top)
    ratioplot.GetLowerRefYaxis().SetTitle(yTitle_bot)
    ratioplot.GetLowerRefYaxis().CenterTitle()
    
    ratioplot.GetLowerRefGraph().SetMinimum(0.01)
    ratioplot.GetLowerRefGraph().SetMaximum(2.00)

    ratioplot.SetGridlines(array.array('d',[1.]),1)
    
    ratioplot.GetLowYaxis().SetNdivisions(505)
    #ratioplot.GetLowYaxis().SetTickLength(0.04)

    legend.Draw("same")

    if addCMSLabel:
        label_cms, label_preliminary, label_luminosity = getLabels_CMS_preliminary(0.1, 0.95, 0.58)
        label_cms.Draw()
        label_preliminary.Draw()
        label_luminosity.Draw()
    
    ratioplot.GetUpperPad().cd()

    hs.Draw("histsame")
    
    h_unc = h_pred.Clone("h_uncertainty")
    setUncertaintyStyle(h_unc)
    legend.AddEntry(h_unc, "Uncertainty", "f")
    h_unc.Draw("e2same")

    h_ratiounc = h_unc.Clone("h_ratiouncertainty")
    
    for i in range(1, h_unc.GetNbinsX()+1):
        bincontent=h_unc.GetBinContent(i)
        binerror=h_unc.GetBinError(i)
        h_ratiounc.SetBinContent(i, 1.)
        if bincontent > 0.:
            h_ratiounc.SetBinError(i, binerror/bincontent)
        else:
            h_ratiounc.SetBinError(i, 0.)
            
    setUncertaintyStyle(h_ratiounc)

    ratioplot.GetLowerPad().cd()

    h_ratiounc.Draw("e2same")
    
    canvas.SaveAs(outname)
            
        
        
def DrawStackHistograms(histograms, outname, blind=True, verbose=False):

    #hists_sorted = sorted(histograms, key=lambda h: h.Integral())
    
    hs = THStack("h_stack","")
    h_obs = None
    
    l = TLegend(0.77,0.60,0.90,0.87)

    for hist in histograms:

        hname=hist.GetName()

        if "data_obs" in hname:
            if not blind:
                h_obs = hist
                SetHistogramStyle(h_obs, "Observed")
                l.AddEntry(h_obs, "Observed", "p")
        else:
            SetHistogramStyle(hist)
            label = GetLabelName(hist)
            hs.Add(hist)
            l.AddEntry(hist,label,"f")
            
    h_est = hs.GetStack().Last().Clone()

    if verbose:
        est_err = Double(0.)
        nbins_est = h_est.GetNbinsX()
        est_yields = h_est.IntegralAndError(1,nbins_est,est_err)
        print "predicted yields :", est_yields,"+/-",est_err
        if not blind:
            print "observed yields :", h_obs.Integral()

    # draw
    
    # SetCMSStyle()   TODO

    canvas = TCanvas()
    gPad.SetLogy()
    
    #ymax_hs = hs.GetMaximum()
    #ymax = ymax_hs if h_obs is None else max(ymax_hs, h_obs.GetMaximum())
    ymax=10000
    ymin=0.1
    
    hs.SetMaximum(ymax)
    hs.SetMinimum(ymin)
    hs.Draw("HIST")
    l.Draw("same")
    if not blind:
        h_obs.Draw("same")

    canvas.SaveAs(outname)

    
            
def SetHistogramStyle(h, label=None):

    if label is None:
        label = h.GetName()
        
    if 'obs' in label.lower():
        h.SetMarkerStyle(20)
        h.SetMarkerColor(1)
    elif 'tth' in label.lower():
        h.SetFillColor(kRed)
        h.SetLineColor(kBlack)
    elif 'ttw' in label.lower():
        h.SetFillColor(kGreen+2)
        h.SetLineColor(kBlack)
    elif 'ttz' in label.lower():
        h.SetFillColor(kGreen-8)
        h.SetLineColor(kBlack)
    elif 'ewk' in label.lower() or 'electroweak' in label.lower():
        h.SetFillColor(kViolet-2)
        h.SetLineColor(kBlack)
    elif 'rare' in label.lower():
        h.SetFillColor(kAzure+1)
        h.SetLineColor(kBlack)
    elif 'fake' in label.lower():
        h.SetFillColor(kBlack)
        h.SetLineColor(kBlack)
        h.SetFillStyle(3335)
    #elif 'flip' in label.lower():
    #    
    else:
        print "WARNING: undefined label"
    

def GetLabelName(hist):
    hname = hist.GetName().lower()

    if 'data_obs' in hname:
        return "Observed"
    elif 'tth' in hname:
        return "ttH"
    elif 'ttz' in hname:
        return "ttZ"
    elif 'ttw' in hname:
        return "ttW"
    elif "ewk" in hname:
        return "Electroweak"
    elif "rares" in hname:
        return "Rares"
    elif "fakes" in hname:
        return "Fakes"
    else:
        return hist.GetName()


def DrawHistfromTree(tree, variable, selection=""):
    tree.Draw(variable+">>htree", selection)
    htree = gDirectory.Get("htree")
    return htree

#def setTDRStyle():
