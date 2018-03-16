from ROOT import THStack, TLegend, TH1, TCanvas, gPad
from ROOT import kRed, kBlack, kGreen, kViolet, kAzure
from ROOT import Double


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
