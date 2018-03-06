import argparse
from ROOT import TCanvas, TFile, TTree, TH1F, TH2F
from ROOT import gDirectory, gStyle

# tau decayMode enum
# 0: kElectron; 1: kMuon; 2: kOneProng0pi0; 3: kOneProng1pi0; 4: kOneProng2pi0;
# 5: kThreeProng0pi0; 6: kThreeProng1pi0; 7: kOther
decayMode_dict = {0:'electron', 1:'muon', 2:'1prong0pi0', 3:'1prong1pi0',
                  4:'1prong2pi0', 5:'3prong0pi0', 6:'3prong1pi0', 7:'other'}

variablelist = ['cos','evis','evisfrac','upsilon','energyasym','cosTauTheta',
                'x1ThreeProngs','x2ThreeProngs','x3ThreeProngs','cosPsi','vismass',
                'evisfrac_mctau','xAsym','evisdiff','evissum']
charges = ['plus','minus']
signs={'plus':'+', 'minus':'-'}
frames = ['lab','bRF','visRF']
framenames = {'lab':'lab', 'bRF':'boson', 'visRF':'visible'}


def getAxisTitleName(variable, charge):
    if variable == 'cos':
        return "cos#theta^{"+signs[charge]+"}"
    elif variable == 'evis':
        return "2E_{vis}^{"+signs[charge]+"}/M_{x}"
    elif variable == 'evisfrac':
        return "E_{vis}^{#tau^{"+signs[charge]+"}}/(E_{vis})"
    #    return "E_{vis}/E_{#tau^{"+signs[charge]+"}}"
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
        return "M_{vis}"
    elif variable == 'evisAsym':
        return "|E_{vis}^{#tau+}-E_{vis}^{#tau-}|/(E_{vis}^{#tau+}+E_{vis}^{#tau-})"
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
    elif variable == 'xAsym':
        x1='x1ThreeProngs_'+charge+'_'+frame
        x2='x2ThreeProngs_'+charge+'_'+frame
        x3='x3ThreeProngs_'+charge+'_'+frame
        name = "(2*min("+x1+","+x2+")-"+x3+")/(2*min("+x1+","+x2+")+"+x3+")"
        #name = "(2*min("+x1+","+x2+")-"+x3+")"
        #print name
    elif variable in ['evisdiff','evissum']:
        name = variable+'_'+frame
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

def makeVisMassHistfromTree(tree, decaymodes, nbins=20, minpt=None, maxeta=None):
    bname = "vismass"
    cuts = getCutsString(decaymodes,'minus',minpt,maxeta)+"&&"
    cuts += getCutsString(decaymodes,'plus',minpt,maxeta)

    h = getHistfromTree(tree, bname, cuts)
    # rebin
    rebinHistogram1D(h, nbins)

    # titles
    h.SetTitle(getDecayModesName(decaymodes,'+'))
    h.GetXaxis().SetTitle(getAxisTitleName('vismass',''))

    h.SetName('vismass_'+getDecayModesName(decaymodes,''))

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
        title += " ("+framenames[frame]+" frame)"
    h2d.SetTitle(title)
    h2d.GetXaxis().SetTitle(getAxisTitleName(variable, 'minus'))
    h2d.GetYaxis().SetTitle(getAxisTitleName(variable, 'plus'))

    # histogram name
    hname = variable+'_corr_'+getDecayModesName(decaymodesX,'')+'_'+getDecayModesName(decaymodesY,'')
    if frame is not None:
        hname += '_'+frame
    h2d.SetName(hname)

    return h2d


def makeEvisAsymHistfromTree(tree, decaymodesX, decaymodesY, frame, nbins=20,
                             selection="", minpt=None, maxeta=None):
    #bdiff = getBranchName('evisdiff', '', frame)
    #bsum = getBranchName('evissum', '', frame)
    #bname = "fabs("+bdiff+")/"+bsum
    
    ###################
    bname1="evis_plus_"+frame
    bname2="evis_minus_"+frame
    bname = "fabs(("+bname1+"-"+bname2+"))/("+bname1+"+"+bname2+")"
    #bname = "(2*"+bname1+"-1)*(2*"+bname2+"-1)"
    ##################
    
    # cuts
    cutsX = getCutsString(decaymodesX, 'minus', minpt, maxeta, selection)
    cutsY = getCutsString(decaymodesY, 'plus', minpt, maxeta, selection)
    cuts = cutsX+"&&"+cutsY
    
    # get histogram
    h = getHistfromTree(tree, bname, cuts)

    # rebin
    rebinHistogram1D(h, nbins)

    # histogram title and axis title
    # Y vs X
    title = getDecayModesName(decaymodesY)+' vs '+getDecayModesName(decaymodesX)
    if frame is not None:
        title += " ("+framenames[frame]+" frame)"
    h.SetTitle(title)
    h.GetXaxis().SetTitle(getAxisTitleName('evisAsym',''))

    # histogram name
    hname = 'evisAsym_'+getDecayModesName(decaymodesX,'')+'_'+getDecayModesName(decaymodesY,'')
    if frame is not None:
        hname += '_'+frame
    h.SetName(hname)

    return h
    

def makeDecayModes2DHistfromTree(tree, decaymodes=[0,1,2,3,5]):
    # cuts
    cuts=getCutsString(decaymodes,'plus',minpt=20.,maxeta=2.4)
    cuts+="&&"+getCutsString(decaymodes,'minus',minpt=20.,maxeta=2.4)
    
    bname = "decayMode_plus:decayMode_minus"
    h2d = getHistfromTree(tree, bname, cuts)

    # normalize
    total=h2d.Integral()
    print total
    h2d.Scale(1./total);
    
    # titles
    h2d.SetTitle("#tau decay modes")
    h2d.GetXaxis().SetTitle("#tau^{-} decay mode")
    h2d.GetYaxis().SetTitle("#tau^{+} decay mode")

    # histogram name
    h2d.SetName("tauDecayModes")
    if 0 not in decaymodes and 1 not in decaymodes:
        h2d.SetName("tauDecayModes_hadronic")

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

def drawMultiHists(hlist, dir='./', optionlist=None):
    assert(len(hlist)>0)
    gStyle.SetOptStat(10)
    canvas = TCanvas()
    hname = hlist[0].GetName()

    for i, h in enumerate(hlist):
        if optionlist is None:
            if i>0:
                h.Draw("same")
            else:
                h.Draw()
        else:
            assert(len(optionlist)==len(hlist))
            h.Draw(optionlist[i])

    canvas.SaveAs(dir+hname+'.pdf')
    
###############################################
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

    decayModesComb = [([2],[2]),([2],[3]),([2],[5]),([3],[3]),([3],[5]),([5],[5]),
                      ([2,3,5],[2,3,5])]

    ###########
    # decaymode
    gStyle.SetPaintTextFormat("1.2f")
    htmp=makeDecayModes2DHistfromTree(tree,decaymodes=[0,1,2,3,5])
    htmp2=makeDecayModes2DHistfromTree(tree,decaymodes=[0,1,2,3,5])
    drawMultiHists([htmp,htmp2],dir=args.outdir,optionlist=["colz","text same"])
    #drawHistogram(htmp, dir=args.outdir,options='colz')
    # hadronic
    htmp=makeDecayModes2DHistfromTree(tree, decaymodes=[2,3,5])
    htmp2=makeDecayModes2DHistfromTree(tree, decaymodes=[2,3,5])
    #drawHistogram(htmp, dir=args.outdir,options='colz')
    drawMultiHists([htmp,htmp2],dir=args.outdir,optionlist=["colz","text same"])

    ###########
    # vis mass
    htmp=makeVisMassHistfromTree(tree, decaymodes=[2,3,5], minpt=20., maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir)
    
    ###########
    # evis aysmmetry
    # 1D
    for (modeX, modeY) in decayModesComb:
        htmp=makeEvisAsymHistfromTree(tree,modeX,modeY,'lab',minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
        htmp=makeEvisAsymHistfromTree(tree,modeX,modeY,'bRF',minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
        htmp=makeEvisAsymHistfromTree(tree,modeX,modeY,'visRF',minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)
    
    ###########
    # evisfrac
    # 1D
    for decays in [[2],[3],[5],[2,3,5]]:
        for charge in ['minus','plus']:
            htmp = makeHist1DfromTree(tree,'evisfrac',charge,decays,'lab',
                                      minpt=20.,maxeta=2.4)
            drawHistogram(htmp, dir=args.outdir)

    for decays in [[2],[3],[5],[2,3,5]]:
        for charge in ['minus','plus']:
            htmp = makeHist1DfromTree(tree,'evisfrac',charge,decays,'bRF',
                                      minpt=20.,maxeta=2.4)
            drawHistogram(htmp, dir=args.outdir)

    #for decays in [[2],[3],[5],[2,3,5]]:
    #    for charge in ['minus','plus']:
    #        htmp = makeHist1DfromTree(tree,'evisfrac_mctau',charge,decays,'lab',
    #                                  minpt=20.,maxeta=2.4)
    #        drawHistogram(htmp, dir=args.outdir)

    #for decays in [[2],[3],[5],[2,3,5]]:
    #    for charge in ['minus','plus']:
    #        htmp = makeHist1DfromTree(tree,'evisfrac_mctau',charge,decays,'bRF',
    #                                  minpt=20.,maxeta=2.4)
    #        drawHistogram(htmp, dir=args.outdir)           

    # 2D
    for (modeX, modeY) in decayModesComb:
        htmp=makeHist2DfromTree(tree,'evisfrac',modeX,modeY,'lab',
                                minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir, options='colz')
        htmp=makeHist2DfromTree(tree,'evisfrac',modeX,modeY,'bRF',
                                minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir, options='colz')
        
        #htmp=makeHist2DfromTree(tree,'evisfrac_mctau',modeX,modeY,'lab',minpt=20.,maxeta=2.4)
        #drawHistogram(htmp, dir=args.outdir, options='colz')
        #htmp=makeHist2DfromTree(tree,'evisfrac_mctau',modeX,modeY,'bRF',minpt=20.,maxeta=2.4)
        #drawHistogram(htmp, dir=args.outdir, options='colz')      

    ###########
    # upsilon
    # 1D
    for decays in [[2],[3],[5],[2,3,5]]:
        for charge in ['minus','plus']:
            htmp = makeHist1DfromTree(tree,'upsilon',charge,decays,'lab',
                                      minpt=20.,maxeta=2.4)
            drawHistogram(htmp, dir=args.outdir)

    # 2D
    for (modeX, modeY) in decayModesComb:
        htmp=makeHist2DfromTree(tree,'upsilon',modeX,modeY,'lab',
                                minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir, options='colz')

    ###########
    # energyasym
    # 1D
    for charge in ['minus','plus']:  # kOneProng1pi0 only
        htmp = makeHist1DfromTree(tree,'energyasym',charge,[3],'lab',
                                  minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)

    # 2D
    for (modeX, modeY) in decayModesComb:
        htmp=makeHist2DfromTree(tree,'energyasym',modeX,modeY,'lab',
                                minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir, options='colz')

    ###########
    # x1ThreeProngs, x2ThreeProngs, x3ThreeProngs
    # 1D
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

    # 2D
    # kThreeProng0pi0 only
    htmp=makeHist2DfromTree(tree,'x3ThreeProngs',[5],[5],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')

    ########### 
    # cosPsi
    # 1D
    for charge in ['minus','plus']:  # kThreeProng0pi0 only
        htmp = makeHist1DfromTree(tree,'cosPsi',charge,[5],minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)

    # 2D
    # kThreeProng0pi0
    htmp=makeHist2DfromTree(tree,'cosPsi',[5],[5],minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
    
    ###########    
    # xAsym: (2*min(x1,x2)-x3) / (2*min(x1,x2)+x3)
    # 1D
    for charge in ['minus','plus']:  # kThreeProng0pi0 only
        htmp = makeHist1DfromTree(tree,'xAsym',charge,[5],'lab',minpt=20.,maxeta=2.4)
        drawHistogram(htmp, dir=args.outdir)

    # 2D
    htmp=makeHist2DfromTree(tree,'xAsym',[5],[5],'lab',minpt=20.,maxeta=2.4)
    drawHistogram(htmp, dir=args.outdir, options='colz')
