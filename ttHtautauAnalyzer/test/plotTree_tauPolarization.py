import argparse
from ROOT import TCanvas, TFile, TTree, TH1F, TH2F
from ROOT import gDirectory, gStyle

# tau decayMode enum
# 0: kElectron; 1: kMuon; 2: kOneProng0pi0; 3: kOneProng1pi0; 4: kOneProng2pi0;
# 5: kThreeProng0pi0; 6: kThreeProng1pi0; 7: kOther
decayMode_dict = {0:'electron', 1:'muon', 2:'1prong0pi0', 3:'1prong1pi0',
                  4:'1prong2pi0', 5:'3prong0pi0', 6:'3prong1pi0', 7:'other'}

variablelist = ['cos','evis','evisfrac','upsilon','energyasym','cosTauTheta',
                'x1ThreeProngs','x2ThreeProngs','x3ThreeProngs','cosPsi','vismass']
charges = ['plus','minus']
signs={'plus':'+', 'minus':'-'}
frames = ['lab','bRF','vis']
framenames = {'lab':'lab', 'bRF':'boson', 'vis':'visible'}


def getAxisTitleName(variable, charge):
    if variable == 'cos':
        return "cos#theta^{"+signs[charge]+"}"
    elif variable == 'evis':
        return "2E_{vis}^{"+signs[charge]+"}/M_{x}"
    elif variable == 'evisfrac':
        return "E_{vis}/E_{#tau^{"+signs[charge]+"}}"
    elif variable == 'upsilon':
        return "2E_{ldgtrk}/E_{vis}-1 (#tau^{"+signs[charge]+"})"
    elif variable == 'energyasym':
        return "(E_{#pm}-E_{0})/(E_{#pm}+E_{0}) (#tau^{"+signs[charge]+"})"
    elif variable == 'cosTauTheta':
        return "Cosine of the angle between #tau^{"+signs[charge]+"} and its visible decay products"
    elif variable in ['x1ThreeProngs','x2ThreeProngs','x3ThreeProngs']:
        return variable[0:2]
    elif variable == 'cosPsi':
        return "cos#Psi_{#tau^{"+signs[charge]+"}}"
    elif variable == 'vismass':
        return "visible mass"
    else:
        return variable+'_'+charge

    
def getDecayModesName(decaymodes, delimiter='+'):
    # special cases
    if sorted(decaymodes)==[0,1]:
        return "leptonic"
    if sorted(decaymodes)==[2,3,5]:
        return "recoOldDM"
    #if sorted(decaymodes)==[2,3,4,5]:  #?
    #    return "recoOldDM"
    
    name=""
    for i, d in enumerate(decaymodes):
        assert(d>=0 and d<=7)
        if i>0:
            name += delimiter
        name += decayMode_dict[d]
        
    return name


def getBranchName(variable, charge, frame):
    name = ""
    if variable == 'vismass':
        name = variable
    elif variable in ['cos','cosPsi']:
        assert(charge is not None)
        name = variable+'_'+charge
    else:
        assert(charge is not None and frame is not None)
        name = variable+'_'+charge+'_'+frame

    return name


def getCutsString(decaymodes, charge, minpt, maxeta, extra=""):
    cuts = extra  # additional cuts
    
    # pt and eta cuts
    if minpt is not None:
        cuts += "&&pt_vis_"+charge+">"+str(minpt)
    if maxeta is not None:
        cuts += "&&fabs(eta_vis_"+charge+")<"+str(maxeta)
        
    #decay modes
    modes="("
    for mode in decaymodes:
        assert(mode>=0 and mode<=7)
        modes += "decayMode_"+charge+"=="+str(mode)
        if mode != decaymodes[-1]:
            modes+="||"
    modes += ')'
    cuts += '&&'+modes
    
    return cuts.lstrip("&")


def getHistfromTree(tree, variables, selections):
    tree.Draw(variables+">>htree", selections)
    htree = gDirectory.Get("htree")
    return htree


def rebinHistogram1D(hist, nbins, verbose=False):
    #hist.Sumw2()
    nbins_old = hist.GetNbinsX()
    if verbose:
        print "current number of bins :", nbins_old
        print "target number of bins :", nbins
        print "ngroup :", nbins_old/nbins
    hist.Rebin(nbins_old/nbins)


def rebinHistogram2D(hist, nbinsx, nbinsy):
    #hist.Sumw2()
    nbinsx_old = hist.GetNbinsX()
    nbinsy_old = hist.GetNbinsY()
    hist.Rebin2D(nbinsx_old/nbinsx, nbinsy_old/nbinsy)
    
    
def makeHist1DfromTree(tree, variable, charge, decaymodes, frame=None, nbins=20,
                           selections="", minpt=None, maxeta=None):
    assert(variable in variablelist)
    assert(charge in charges or charge is None)
    assert(frame in frames or frame is None)
    
    # branch name in tree
    bname = getBranchName(variable, charge, frame)

    # cuts
    cuts = getCutsString(decaymodes, charge, minpt, maxeta, selections)

    # get histogram
    h = getHistfromTree(tree, bname, cuts)

    # rebin
    rebinHistogram1D(h, nbins)
    
    # histogram title and axis title
    title = "#tau^{"+signs[charge]+"}  "      
    title += getDecayModesName(decaymodes,'+')
    if frame is not None:
        title += "  "+framenames[frame]+' frame'
    
    h.SetTitle(title)
    h.GetXaxis().SetTitle(getAxisTitleName(variable, charge))

    # histogram name
    hname = variable+'_'+charge+'_'+getDecayModesName(decaymodes,'')
    if frame is not None:
        hname += '_'+frame 
    h.SetName(hname)
    
    return h


def makeHist2DfromTree(tree, variable, decaymodesX, decaymodesY, frame=None,
                       nbinsx=20, nbinsy=20, selection="", minpt=None, maxeta=None):
    assert(variable in variablelist)  # TODO: variable based on decay mode
    assert(variable != 'vismass')
    assert(frame in frames or frame is None)

    # branch names in tree
    # tau minus on x-axis, tau plus on y-axis
    bnameX = getBranchName(variable, 'minus', frame)
    bnameY = getBranchName(variable, 'plus', frame)
    bname = bnameY+':'+bnameX

    # cuts
    cutsX = getCutsString(decaymodesX, 'minus', minpt, maxeta, selection)
    cutsY = getCutsString(decaymodesY, 'plus', minpt, maxeta, selection)
    cuts = cutsX+"&&"+cutsY

    # get histogram
    h2d = getHistfromTree(tree, bname, cuts)

    # rebin
    rebinHistogram2D(h2d, nbinsx, nbinsy)

    # histogram title and axis title
    # Y vs X
    title = getDecayModesName(decaymodesY)+' vs '+getDecayModesName(decaymodesX)
    if frame is not None:
        title += "("+framenames[frame]+" frame)"
    h2d.SetTitle(title)
    h2d.GetXaxis().SetTitle(getAxisTitleName(variable, 'minus'))
    h2d.GetYaxis().SetTitle(getAxisTitleName(variable, 'plus'))

    # histogram name
    hname = variable+'_corr_'+getDecayModesName(decaymodesX,'')+'_'+getDecayModesName(decaymodesY,'')
    if frame is not None:
        hname += '_'+frame
    h2d.SetName(hname)

    return h2d
    

def drawHistogram(h, dir='./', options=None):
    gStyle.SetOptStat(10)
    canvas = TCanvas()
    hname = h.GetName()
    
    if options is None:
        h.Draw()
    else:
        h.Draw(options)
        
    canvas.SaveAs(dir+hname+'.pdf')
    #canvas.SaveAs(dir+hname+'.png')
    

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help="input root file name")
    parser.add_argument('-t','--treename', type=str, default='GenAna/eventTree',
                        help="tree name")
    parser.add_argument('-o',"--outdir",type=str, default='ttH/', help="output directory")
    
    args = parser.parse_args()

    # read root file and get event tree
    fin = TFile(args.infile, 'read')
    tree = fin.Get(args.treename)
    
    # tau decayMode enum
    # 0: kElectron; 1: kMuon; 2: kOneProng0pi0; 3: kOneProng1pi0; 4: kOneProng2pi0;
    # 5: kThreeProng0pi0; 6: kThreeProng1pi0; 7: kOther

    ### 1D
    # evisfrac
    for decays in [[2],[3],[5],[2,3,5]]:
        for charge in ['minus','plus']:
            htmp = makeHist1DfromTree(tree,'evisfrac',charge,decays,'lab',
                                      minpt=20.,maxeta=2.4)
            drawHistogram(htmp, dir=args.outdir)
    
    # upsilon
    for decays in [[2],[3],[5],[2,3,5]]:
        for charge in ['minus','plus']:
            htmp = makeHist1DfromTree(tree,'upsilon',charge,decays,'lab',
                                      minpt=20.,maxeta=2.4)
            drawHistogram(htmp, dir=args.outdir)
    # energyasym
    for charge in ['minus','plus']:  # kOneProng1pi0 only
        htmp = makeHist1DfromTree(tree,'energyasym',charge,[3],'lab',
                                  minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
        
    # x1ThreeProngs, x2ThreeProngs, x3ThreeProngs
    for charge in ['minus','plus']:  # kThreeProng0pi0 only
        htmp = makeHist1DfromTree(tree,'x1ThreeProngs',charge,[5],'lab',
                                  minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
        htmp = makeHist1DfromTree(tree,'x2ThreeProngs',charge,[5],'lab',
                                  minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
        htmp = makeHist1DfromTree(tree,'x3ThreeProngs',charge,[5],'lab',
                                  minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)

    # cosPsi
    for charge in ['minus','plus']:  # kThreeProng0pi0 only
        htmp = makeHist1DfromTree(tree,'cosPsi',charge,[5],minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)

    # vismass
    # TODO
        
    ### 2D
    # evisfrac
    # kOneProng0pi0
    htmp=makeHist2DfromTree(tree,'evisfrac',[2],[2],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
    # kOneProng1pi0
    htmp=makeHist2DfromTree(tree,'evisfrac',[3],[3],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
    # kThreeProng0pi0
    htmp=makeHist2DfromTree(tree,'evisfrac',[5],[5],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
        
    # upsilon
    # kOneProng0pi0
    htmp=makeHist2DfromTree(tree,'upsilon',[2],[2],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
    # kOneProng1pi0
    htmp=makeHist2DfromTree(tree,'upsilon',[3],[3],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
    # kThreeProng0pi0
    htmp=makeHist2DfromTree(tree,'upsilon',[5],[5],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')

    # energyasym
    # kOneProng1pi0
    htmp=makeHist2DfromTree(tree,'energyasym',[3],[3],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')

    # cosPsi
    # kThreeProng0pi0
    htmp=makeHist2DfromTree(tree,'cosPsi',[5],[5],minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
