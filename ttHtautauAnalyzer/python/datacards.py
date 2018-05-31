from ROOT import TTree, TH1D, TH2D
from math import sqrt
from plotting import DrawHistfromTree

DataSamples = ['DoubleMuon','SingleMuon','DoubleEG','SingleElectron','MuonEG']

SamplesInChannel = {'ttH':['ttH'],
                    'TTW':['TTW'],
                    'TTZ':['TTZ','TTGJets'],
                    'EWK':['WZ','ZZ','WW','WG','ZG'],
                    'Rares':['TTTT','tZq','WWW','WWZ','WZZ','ZZZ','WpWp'],
                    'fakes_data':DataSamples,
                    'flips_data':DataSamples,
                    'data_obs':DataSamples}

datasamples = ['data_e','data_mu','data_dieg','data_dimu','data_mueg']
SamplesInChannel2017 = {'ttH':['ttHJetToNonbb'],
                        'TTW':['TTW_psw','TTWW'],
                        'TTZ':['TTZ','TTZ_M1to10'],
                        'EWK':['WZTo3LNu','ZZ_ext','WW'],
                        'Rares':['WWW','WWZ','WZZ','ZZZ','tZq','TTTT','WpWpJJ', 'WWTo2L2Nuds'],
                        'tH':['tHq','tHW'],
                        'Conversion':['TTGJets'],
                        'ggH':['ggHZZ4l'],
                        'VH':['VHToNonbb'],
                        'TT':['TTToDiLep','TTToSemiLep','TTToHad'],
                        'ST':['ST_sLep_psw','ST_tT','ST_tTbar','ST_tWT_psw','ST_tWTbar_psw','ST_tWll'],
                        'DY':['DYJets_M50','DYJets_M50_ext','DYJets_M10to50'],
                        'fakes_data':datasamples,
                        'flips_data':datasamples,
                        'data_obs':datasamples} 

BTagSysts = ['LFUp','LFDown','HFUp','HFDown',
             'HFStats1Up','HFStats1Down','HFStats2Up','HFStats2Down',
             'LFStats1Up','LFStats1Down','LFStats2Up','LFStats2Down',
             'cErr1Up','cErr1Down','cErr2Up','cErr2Down']

ThSysts = ['x1Up','x1Down','y1Up','y1Down']

FakeRateLepSysts = ['FRe_normUp','FRe_normDown','FRe_ptUp','FRe_ptDown',
                    #'FRe_bUp','FRe_bDown','FRe_ecUp','FRe_ecDown',
                    'FRe_beUp','FRe_beDown',
                    'FRm_normUp','FRm_normDown','FRm_ptUp','FRm_ptDown',
                    #'FRm_bUp','FRm_bDown','FRm_ecUp','FRm_ecDown',
                    'FRm_beUp','FRm_beDown']

FakeTauSysts = ['FRjt_normUp','FRjt_normDown','FRjt_shapeUp','FRjt_shapeDown']

ClosureTests = ['Clos_e_shapeUp','Clos_e_shapeDown','Clos_m_shapeUp',
                'Clos_m_shapeDown']

def getEventWeight(event, name):

    weight = event.event_weight

    if 'btag_LFUp' in name:
        weight = event.event_weight_btag_LFUp
    elif 'btag_LFDown' in name:
        weight = event.event_weight_btag_LFDown
    elif 'btag_HFUp' in name:
        weight = event.event_weight_btag_HFUp
    elif 'btag_HFDown' in name:
        weight = event.event_weight_btag_HFDown
    elif 'btag_HFStats1Up' in name:
        weight = event.event_weight_btag_HFStats1Up
    elif 'btag_HFStats1Down' in name:
        weight = event.event_weight_btag_HFStats1Down
    elif 'btag_HFStats2Up' in name:
        weight = event.event_weight_btag_HFStats2Up
    elif 'btag_HFStats2Down' in name:
        weight = event.event_weight_btag_HFStats2Down
    elif 'btag_LFStats1Up' in name:
        weight = event.event_weight_btag_LFStats1Up
    elif 'btag_LFStats1Down' in name:
        weight = event.event_weight_btag_LFStats1Down
    elif 'btag_LFStats2Up' in name:
        weight = event.event_weight_btag_LFStats2Up
    elif 'btag_LFStats2Down' in name:
        weight = event.event_weight_btag_LFStats2Down
    elif 'btag_cErr1Up' in name:
        weight = event.event_weight_btag_cErr1Up
    elif 'btag_cErr1Down' in name:
        weight = event.event_weight_btag_cErr1Down
    elif 'btag_cErr2Up' in name:
        weight = event.event_weight_btag_cErr2Up
    elif 'btag_cErr2Down' in name:
        weight = event.event_weight_btag_cErr2Down
    elif 'thu_shape_x1Up' in name:
        weight = event.event_weight_thu_shape_x1Up
    elif 'thu_shape_x1Down' in name:
        weight = event.event_weight_thu_shape_x1Down
    elif 'thu_shape_y1Up' in name:
        weight = event.event_weight_thu_shape_y1Up
    elif 'thu_shape_y1Down' in name:
        weight = event.event_weight_thu_shape_y1Down
    elif 'FRjt_normUp' in name:
        weight = event.event_weight_FRjt_normUp
    elif 'FRjt_normDown' in name:
        weight = event.event_weight_FRjt_normDown
    elif 'FRjt_shapeUp' in name:
        weight = event.event_weight_FRjt_shapeUp
    elif 'FRjt_shapeDown' in name:
        weight = event.event_weight_FRjt_shapeDown
    elif 'FRe_normUp' in name:
        weight = event.event_weight_FRe_normUp
    elif 'FRe_normDown' in name:
        weight = event.event_weight_FRe_normDown
    elif 'FRe_ptUp' in name:
        weight = event.event_weight_FRe_ptUp
    elif 'FRe_ptDown' in name:
        weight = event.event_weight_FRe_ptDown
    elif 'FRe_bUp' in name:
        weight = event.event_weight_FRe_bUp
    elif 'FRe_bDown' in name:
        weight = event.event_weight_FRe_bDown
    elif 'FRe_ecUp' in name:
        weight = event.event_weight_FRe_ecUp
    elif 'FRe_ecDown' in name:
        weight = event.event_weight_FRe_ecDown
    elif 'FRm_normUp' in name:
        weight = event.event_weight_FRm_normUp
    elif 'FRm_normDown' in name:
        weight = event.event_weight_FRm_normDown
    elif 'FRm_ptUp' in name:
        weight = event.event_weight_FRm_ptUp
    elif 'FRm_ptDown' in name:
        weight = event.event_weight_FRm_ptDown
    elif 'FRm_bUp' in name:
        weight = event.event_weight_FRm_bUp
    elif 'FRm_bDown' in name:
        weight = event.event_weight_FRm_bDown
    elif 'FRm_ecUp' in name:
        weight = event.event_weight_FRm_ecUp
    elif 'FRm_ecDown' in name:
        weight = event.event_weight_FRm_ecDown

    return weight
        

def getShapeFromTree(tree, histname, var_x, var_y, nbin, xmin, xmax, binningMap):
    
    h = TH1D(histname,"", nbin, xmin, xmax)
    h.Sumw2()
    
    for ev in tree:
        ibin = getBin(ev.GetLeaf(var_x).GetValue(), ev.GetLeaf(var_y).GetValue(),
                      binningMap)

        if 'gentau' in histname.lower() and not ev.isGenMatchedTau:
            continue
        if 'faketau' in histname.lower() and ev.isGenMatchedTau:
            continue
        if 'tth_htt' in histname.lower() and abs(ev.HiggsDecayType)!=15:
            continue
        if 'tth_hww' in histname.lower() and abs(ev.HiggsDecayType)!=24:
            continue
        if 'tth_hzz' in histname.lower() and abs(ev.HiggsDecayType)!=23:
            continue
    
        weight = getEventWeight(ev, histname)
        h.Fill(ibin, weight)
        
    return h

def getShapeFromMergingTrees(trees, histname, var_x, var_y, nbin, xmin, xmax,
                             binningMap):
    # histogram
    h = TH1D(histname,"", nbin, xmin, xmax)
    h.Sumw2()
    eventIDlist = []
    
    for tree in trees:
        for ev in tree:
            eventID = (ev.run, ev.lumi, ev.nEvent)
            #eventID = (ev.run, int(ev.GetLeaf('ls').GetValue()), ev.nEvent)
            
            if eventID in eventIDlist:
                continue
                
            eventIDlist.append(eventID)

            ibin=getBin(ev.GetLeaf(var_x).GetValue(),ev.GetLeaf(var_y).GetValue(),
                        binningMap)
            weight = getEventWeight(ev, histname)

            h.Fill(ibin, weight)

    #print eventIDlist
    return h  # also eventIDlist?

def getBin(x, y, binningMap):
    if binningMap is None:
        return x
    else:
        return read2DHist(x, y, binningMap)
            

def read2DHist(x, y, h2d):
    
    xaxis = h2d.GetXaxis()
    nbinx = xaxis.GetNbins()
    xbin = min(nbinx, max(1,xaxis.FindBin(x)))

    yaxis = h2d.GetYaxis()
    nbiny = yaxis.GetNbins()
    ybin = min(nbiny, max(1,yaxis.FindBin(y)))

    return h2d.GetBinContent(xbin, ybin)

def compIntegral(histogram, includeUnderflowBin, includeOverflowBin):

    sumBinContent = 0.
    numBins = histogram.GetNbinsX()
    firstBin = 0 if includeUnderflowBin else 1
    lastBin = numBins+1 if includeOverflowBin else numBins

    for iBin in range(firstBin, lastBin+1):
        sumBinContent += histogram.GetBinContent(iBin)

    return sumBinContent

def makeBinContentsPositive(histogram, verbosity=False):

    if verbosity:
        print '<make>BinContentsPositive>:'
        print ' integral(', histogram.GetName(),') =', histogram.Integral()
    
    integral_original = max(compIntegral(histogram, True, True), 0.)

    if verbosity:
        print ' integral_original =', integral_original

    numBins = histogram.GetNbinsX()

    for iBin in range(numBins+2):
        binContent_original = histogram.GetBinContent(iBin)
        binError2_original = (histogram.GetBinError(iBin))**2

        if binContent_original < 0.:
            binContent_modified = 0.
            binError2_modified = binError2_original + (binContent_original-binContent_modified)**2
            assert(binError2_modified >= 0.)

            if verbosity:
                print "bin #", iBin, ' (x =', histogram.GetBinCenter(iBin),'): binContent =', binContent_original,' +/-', sqrt(binError2_original), ' --> setting it to binContent =', binContent_modified, ' +/-', sqrt(binError2_modified)

            histogram.SetBinContent(iBin, binContent_modified)
            histogram.SetBinError(iBin, sqrt(binError2_modified))

    integral_modified = max(compIntegral(histogram, True, True), 0.)

    if verbosity:
        print ' integral_modified =', integral_modified

    if integral_modified > 0.:
        sf = integral_original/integral_modified

        if verbosity:
            print '--> scaling histogram by factor =', sf

        histogram.Scale(sf)
    else:
        for iBin in range(numBins+2):
            histogram.SetBinContent(iBin, 0.)

    if verbosity:
        print ' integral(', histogram.GetName(), ') =', histogram.Integral()


def getNtupleFileName_mc(ntuplelist, anatype, sample, correction=None):
    # get ntuple file name from list
    target = '_'+sample
    if correction is not None and correction!='':
        assert(correction.lower()=='jesup' or correction.lower()=='jesdown' or
               correction.lower()=='tesup' or correction.lower()=='tesdown')
        target += '_'+correction.lower()
    target += '_'+anatype

    assert(anatype=='1l2tau' or anatype=='2lss1tau' or anatype=='3l1tau')
        
    with open(ntuplelist) as f:
        for line in f:
            if target in line.strip():
                return line.strip()

        print 'WARNING: CANNOT find ntuple file in', ntuplelist,'(',anatype,sample,correction,')'
        print 'Make sure to hadd and produce it first.'

def getNtupleFileName_data(ntuplelist, anatype, channel, sample):
    # get ntuple file name from list

    assert('data_obs' in channel or 'fakes_data' in channel or 'flips_data' in channel)
    assert('SingleMuon' in sample or 'SingleElectron' in sample or 'DoubleMuon' in sample or 'DoubleEG' in sample or 'MuonEG' in sample)
    
    with open(ntuplelist) as f:
        for line in f:
            if anatype in line.strip() and channel in line.strip() and sample in line.strip():
                return line.strip()

        print 'WARNING: CANNOT find ntuple file in', ntuplelist,'(',anatype,channel,sample,')'
        print 'Make sure to hadd and produce it first.'


def getClosureTestShape(h_nominal, closTest, infile, syst_coname='_CMS_ttHl_',
                        hname_ele = 'x_TT_DL_FR_TT_MC_minus_FR_QCD_MC_ele',
                        hname_mu = 'x_TT_DL_FR_TT_MC_minus_FR_QCD_MC_mu'):

    assert(h_nominal.GetName()=='x_fakes_data')
    assert(closTest in ClosureTests)
    
    # get number of bins and x range from h_nominal
    nbin = h_nominal.GetNbinsX()
    xmin = h_nominal.GetBinLowEdge(1)
    xmax = h_nominal.GetBinLowEdge(nbin)+h_nominal.GetBinWidth(nbin)

    # open file and get histograms
    f_clos = TFile(infile, 'read')

    hname = hname_ele if 'Clos_e_' in closTest else hname_mu
    assert('Clos_e' in closTest or 'Clos_m' in closTest)
    
    h_ttbar_minus_qcd_fr = f_clos.Get(hname)

    h_clos = TH1D("x_fakes_data"+syst_coname+closTest)
    h_clos.Sumw2()

    factor = 1. if 'shapeUp' in closTest else -1.

    h_clos.Add(h_nominal, h_ttbar_minus_qcd_fr, 1., factor)

    return h_clos


def printYields(label, histlist, var='x', genTauSplit=True):
    hname1 = var+'_'+label
    hname2 = var+'_'+label+'_gentau'
    hname3 = var+'_'+label+'_faketau'

    yields = [-1., -1., -1.]

    for hist in histlist:
        
        if hname1 == hist.GetName():
            yields[0]=hist.Integral()
        elif hname2 == hist.GetName():
            yields[1]=hist.Integral()
        elif hname3 == hist.GetName():
            yields[2]=hist.Integral()

        if yields[0]!= -1. and yields[1]!=-1. and yields[2]!=-1.:
            break

    print label,'\t','%.5f'%yields[0],
    if genTauSplit:
        print '\t','%.5f'%yields[1],'\t','%.5f'%yields[2]
    else:
        print
