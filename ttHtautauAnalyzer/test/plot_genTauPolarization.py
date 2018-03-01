import argparse
from ROOT import TCanvas, TFile, TTree, TH1F, TH2F
from ROOT import gStyle
from math import sin, cos, asin, acos

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str,
                    help="input root file name")
parser.add_argument('-t','--treename', type=str, default='GenAna/eventTree',
                    help="tree name")
parser.add_argument('-o','--outdir', type=str, default='ttH/',
                    help="output directory")
parser.add_argument('-e','--epsilon', type=float, default=0.000001,
                    help="small addition to xmax to include the upper edge in the last bin; set to 0 if wish to exclude upper edge in the last bin")
args = parser.parse_args()

# read root file and get the event tree
fin = TFile(args.infile, 'read')
tree = fin.Get(args.treename)

# tau decayMode enum
# 0: kElectron; 1: kMuon; 2: kOneProng0pi0; 3: kOneProng1pi0; 4: kOneProng2pi0;
# 5: kThreeProng0pi0; 6: kThreeProng1pi0; 7: kOther

###### histograms ######
decaymodes=['leptonic', '1prong0pi0', '1prong1pi0', '1prong2pi0', '3prong0pi0',
            'combDMold']
frames=['lab', 'boson', 'vis']
charges=['plus', 'minus']
signs={'plus':'+','minus':'-'}

### 1D
# h_evis[frame][tau charge][decay mode]
h_evis = [[[] for j in range(len(charges))] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for decay in decaymodes:
            xmax=50. if frame=='lab' else 1.+args.epsilon
            h_tmp = TH1F('h_evis_'+charge+'_'+decay+"_"+frame,
                         '#tau^{'+signs[charge]+'}  '+decay+'  '+frame+' frame',
                         20, 0., xmax)
            h_tmp.GetXaxis().SetTitle("2E_{vis}^{"+signs[charge]+"}/M_{x}")
            h_evis[f][c].append(h_tmp)

# h_energyasym[frame][tau charge][decay mode]
h_energyasym = [[[] for j in range(len(charges))] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for decay in decaymodes:
            h_tmp = TH1F('h_energyasym_'+charge+'_'+decay+"_"+frame,
                         '#tau^{'+signs[charge]+'}  '+decay+'  '+frame+' frame',
                         20, -1., 1.+args.epsilon)
            h_tmp.GetXaxis().SetTitle("(E_{#pm}-E_{0})/(E_{#pm}+E_{0})")
            h_energyasym[f][c].append(h_tmp)
            
# h_cos[tau charge][decay mode]
h_cos = [[] for i in range(len(charges))]

for c, charge in enumerate(charges):
    for decay in decaymodes:
        h_tmp = TH1F('h_cos_'+charge+'_'+decay,'#tau^{'+signs[charge]+'}  '+decay,
                     20, -1., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle("cos#theta^{"+signs[charge]+"}")
        h_cos[c].append(h_tmp)

# h_cosPsi[tau charge]
h_cosPsi_3prongs = [] 

for charge in charges:
    h_tmp = TH1F('h_cosPsi_'+charge+'_'+decay,'#tau^{'+signs[charge]+'} 3prong0pi0',
                 20, -1., 1.+args.epsilon)
    h_tmp.GetXaxis().SetTitle("cos#Psi^{"+signs[charge]+"}")
    h_cosPsi_3prongs.append(h_tmp)
    
### 2D
# nmode: number of combinatinos of tau pair decay mode
nmode = len(decaymodes)*(len(decaymodes)-1)/2 + 1
combModes = []
for i in range(len(decaymodes)-1):
    for j in range(i, len(decaymodes)-1):
        combModes.append(decaymodes[i]+"_"+decaymodes[j])
combModes.append("oldDM_oldDM")
print combModes

# h_costheta_corr[mode]
h_costheta_corr = []
for mode in combModes:
    h_tmp = TH2F('h_costheta_corr_'+mode+'_'+frame, mode+' ('+frame+' frame)',
                 20, -1., 1.+args.epsilon, 20, -1., 1.+args.epsilon)
    h_tmp.GetXaxis().SetTitle("cos#theta^{-}")
    h_tmp.GetYaxis().SetTitle("cos#theta^{+}")
    h_costheta_corr.append(h_tmp)

# h_cosPsi_corr
h_cosPsi_corr = TH2F('h_cosPsi_corr_3prongs', '3prongs_3prongs',
                     20, -1., 1.+args.epsilon, 20, -1., 1.+args.epsilon)
h_cosPsi_corr.GetXaxis().SetTitle("cos#Psi^{-}")
h_cosPsi_corr.GetYaxis().SetTitle("cos#Psi^{+}")

# h_xThreeProngs[frame][charge]
h_xThreeProngs = [[] for i in range(len(frames))]
for f, frame in enumerate(frames):
    for charge in charges:
        h_tmp = TH2F('h_xThreeProngs_tau'+charge+'_'+frame,
                     '3prongs (tau'+charge+', '+frame+' frame)',
                     20, 0., 1.+args.epsilon, 20, 0., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle("x3")
        h_tmp.GetYaxis().SetTitle("min(x1,x2)")
        h_xThreeProngs[f].append(h_tmp)

# h_evis_corr[frame][mode]
h_evis_corr = [[] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for mode in combModes:
        xmax=50. if frame=='lab' else 1.+args.epsilon
        h_tmp = TH2F('h_evis_corr_'+mode+'_'+frame, mode+' ('+frame+' frame)',
                     20, 0., xmax, 20, 0., xmax)
        h_tmp.GetXaxis().SetTitle('2E_{vis}^{-}/M_{X} ('+mode.split('_')[1]+')')
        h_tmp.GetYaxis().SetTitle('2E_{vis}^{+}/M_{X} ('+mode.split('_')[0]+')')
        h_evis_corr[f].append(h_tmp)

# h_energyasym_corr[frame][mode]
h_energyasym_corr = [[] for i in range(len(frames))]
for f, frame in enumerate(frames):
    for mode in combModes:
        h_tmp = TH2F('h_energyasym_corr_'+mode+'_'+frame, mode+' ('+frame+' frame)',
                     20, -1., 1.+args.epsilon, 20, -1., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle('Energyasym #tau^{-} ('+mode.split('_')[1]+')')
        h_tmp.GetYaxis().SetTitle('Energyasym #tau^{+} ('+mode.split('_')[0]+')')
        h_energyasym_corr[f].append(h_tmp)

##################
# tau Helicity Tag
# h_tauHelTag[frame][tau charge][decay mode]
h_tauHelTag = [[[] for j in range(len(charges))] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for decay in decaymodes:
            h_tmp = TH1F('h_tauHelTag_'+charge+'_'+decay+'_'+frame,
                         '#tau^{'+signs[charge]+'}  '+decay+'  '+frame+' frame',
                         20, -1., 1.+args.epsilon)
            h_tmp.GetXaxis().SetTitle("helTag_{#tau^{"+signs[charge]+"}}")
            h_tauHelTag[f][c].append(h_tmp)

# h_tauHelTag_prod[frame][mode]
# h2d_tauHelTag_correlation[frame][mode]
h_tauHelTag_prod = [[] for i in range(len(frames))]
h2d_tauHelTag_correlation = [[] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for mode in combModes:
        h_tmp = TH1F("h_tauHelTag_prod_"+mode+'_'+frame, mode+' ('+frame+' frame)',
                     20, -1., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle("helTag_{#tau^{+}}*helTag_{#tau^{-}}")
        h_tauHelTag_prod[f].append(h_tmp)

        h2d_tmp = TH2F("h2d_tauHelTag_correlation_"+mode+'_'+frame,
                       mode+' ('+frame+' frame)',
                       20, -1., 1.+args.epsilon, 20, -1., 1.+args.epsilon)
        h2d_tmp.GetXaxis().SetTitle("helTag #tau^{-}")
        h2d_tmp.GetYaxis().SetTitle("helTag #tau^{+}")
        h2d_tauHelTag_correlation[f].append(h2d_tmp)

def tauHelValueMapping_m1top1(x): # e.g. energyasym: [-1,1]
    # return 1 if -1 or 1, return -1 if 0
    #htag = 1-2*sin(acos(x))
    htag = 2*abs(x)-1
    return htag

def tauHelValueMapping_0top1(x): # e.g. x3: [0,1]
    # return 1 if 0 or 1, return -1 if 0.5
    xnew = 2.*x-1.
    return tauHelValueMapping_m1top1(xnew)
        
##################

# loop over events in tree
evcnt = 0
for ev in tree:
    evcnt+=1
    if not evcnt%10000:
        print 'processing event', evcnt
    
    # e/mu merged into one leptonic category
    # only consider 1prong0pi0, 1prong1pi0, 1prong2pi0 and 3prong0pi0
    decaymode_p = max(ev.decayMode_plus-1, 0) 
    decaymode_m = max(ev.decayMode_minus-1, 0)
    # 0: leptonic; 1: 1prong0pi0; 2: 1prong1pi0; 3: 1prong2pi0; 4: 3prong0pi0;

    if decaymode_p > 4:
        continue
    if decaymode_m > 4:
        continue
    
    ## lab frame
    # tau+
    h_evis[0][0][decaymode_p].Fill(ev.evis_plus_lab)
    h_energyasym[0][0][decaymode_p].Fill(ev.energyasym_plus_lab)
    # tau-
    h_evis[0][1][decaymode_m].Fill(ev.evis_minus_lab)
    h_energyasym[0][1][decaymode_m].Fill(ev.energyasym_minus_lab)
    
    ## boson rest frame
    # tau+
    h_evis[1][0][decaymode_p].Fill(ev.evis_plus_bRF)
    h_energyasym[1][0][decaymode_p].Fill(ev.energyasym_plus_bRF)
    # tau-
    h_evis[1][1][decaymode_m].Fill(ev.evis_minus_bRF)
    h_energyasym[1][1][decaymode_m].Fill(ev.energyasym_minus_bRF)
    
    ## visible tau pair rest frame
    # tau+
    h_evis[2][0][decaymode_p].Fill(ev.evis_plus_visRF)
    h_energyasym[2][0][decaymode_p].Fill(ev.energyasym_plus_visRF)
    # tau-
    h_evis[2][1][decaymode_m].Fill(ev.evis_minus_visRF)
    h_energyasym[2][1][decaymode_m].Fill(ev.energyasym_minus_visRF)
    
    h_cos[0][decaymode_p].Fill(ev.cos_plus)
    h_cos[1][decaymode_m].Fill(ev.cos_minus)

    if decaymode_p == 4:
        h_cosPsi_3prongs[0].Fill(ev.cosPsi_plus)
        h_xThreeProngs[0][0].Fill(ev.x3ThreeProngs_plus_lab,
                                  min(ev.x1ThreeProngs_plus_lab,
                                      ev.x2ThreeProngs_plus_lab))  # lab frame
        h_xThreeProngs[1][0].Fill(ev.x3ThreeProngs_plus_bRF,
                                  min(ev.x1ThreeProngs_plus_bRF,
                                      ev.x2ThreeProngs_plus_bRF))  # boson rest frame
    if decaymode_m == 4:
        h_cosPsi_3prongs[1].Fill(ev.cosPsi_minus)
        h_xThreeProngs[0][1].Fill(ev.x3ThreeProngs_minus_lab,
                                  min(ev.x1ThreeProngs_minus_lab,
                                      ev.x2ThreeProngs_minus_lab)) # lab frame
        h_xThreeProngs[1][1].Fill(ev.x3ThreeProngs_minus_bRF,
                                  min(ev.x1ThreeProngs_minus_bRF,
                                      ev.x2ThreeProngs_minus_bRF)) # boson rest frame

    #####################
    heltagminus_lab = -233.
    heltagplus_lab = -233.
    heltagminus_bRF = -233.
    heltagplus_bRF = -233.
    
    # tau+
    if decaymode_p==2 or decaymode_p==3:  # 1prong1pi0, 1prong2pi0
        heltagplus_lab = tauHelValueMapping_m1top1(ev.energyasym_plus_lab)
        heltagplus_bRF = tauHelValueMapping_m1top1(ev.energyasym_plus_bRF)
    elif decaymode_p==4: # 3prong0pi0
        heltagplus_lab = tauHelValueMapping_0top1(ev.x3ThreeProngs_plus_lab)
        heltagplus_bRF = tauHelValueMapping_0top1(ev.x3ThreeProngs_plus_bRF)
    #elif decaymode_p==1: # 1prong0pi0
        # TODO 
    # tau-
    if decaymode_m==2 or decaymode_m==3: # 1prong1pi0, 1prong2pi0
        heltagminus_lab = tauHelValueMapping_m1top1(ev.energyasym_minus_lab)
        heltagminus_bRF = tauHelValueMapping_m1top1(ev.energyasym_minus_bRF)
    elif decaymode_m==4: # 3prong0pi0
        heltagminus_lab = tauHelValueMapping_0top1(ev.x3ThreeProngs_minus_lab)
        heltagminus_bRF = tauHelValueMapping_0top1(ev.x3ThreeProngs_minus_bRF)
    #elif decaymode_m==1: # 1prong0pi0
        # TODO  
    
    h_tauHelTag[0][0][decaymode_p].Fill(heltagplus_lab)
    h_tauHelTag[0][1][decaymode_m].Fill(heltagminus_lab)
    h_tauHelTag[1][0][decaymode_p].Fill(heltagplus_bRF)
    h_tauHelTag[1][1][decaymode_m].Fill(heltagminus_bRF)
        
    # Old DecayMode
    index_oldDM = len(decaymodes)-1
    if decaymode_p in [1,2,4] and decaymode_m in [1,2,4]:
        h_evis[0][0][index_oldDM].Fill(ev.evis_plus_lab)
        h_energyasym[0][0][index_oldDM].Fill(ev.energyasym_plus_lab)
        h_evis[1][0][index_oldDM].Fill(ev.evis_plus_bRF)
        h_energyasym[1][0][index_oldDM].Fill(ev.energyasym_plus_bRF)
        h_evis[2][0][index_oldDM].Fill(ev.evis_plus_visRF)
        h_energyasym[2][0][index_oldDM].Fill(ev.energyasym_plus_visRF)
        
        h_evis[0][1][index_oldDM].Fill(ev.evis_plus_lab)
        h_energyasym[0][1][index_oldDM].Fill(ev.energyasym_plus_lab)
        h_evis[1][1][index_oldDM].Fill(ev.evis_minus_bRF)
        h_energyasym[1][1][index_oldDM].Fill(ev.energyasym_minus_bRF)
        h_evis[2][1][index_oldDM].Fill(ev.evis_minus_visRF)
        h_energyasym[2][1][index_oldDM].Fill(ev.energyasym_minus_visRF)

        h_tauHelTag[0][0][index_oldDM].Fill(heltagplus_lab)
        h_tauHelTag[0][1][index_oldDM].Fill(heltagminus_lab)
        h_tauHelTag[1][0][index_oldDM].Fill(heltagplus_bRF)
        h_tauHelTag[1][1][index_oldDM].Fill(heltagminus_bRF)

    ### Correlations    
    if decaymode_p > decaymode_m:
        continue  # only plot e.g. l+pi-, but not l-pi+
    
    # Index of tau pair decay mode in combModes list
    #(2*len(combModes)-decaymode_p+1)*decaymode_p/2 + decaymode_m - decaymode_p
    imode = (2*(len(decaymodes)-1)-decaymode_p-1) * decaymode_p / 2 + decaymode_m
    #print decaymodes[decaymode_p], decaymodes[decaymode_m]
    #print combModes[imode]
    assert(decaymodes[decaymode_p]+'_'+decaymodes[decaymode_m]==combModes[imode])

    ## lab frame
    #varhelicity_plus_lab = ev.evis_plus_lab
    #if decaymodes[decaymode_p]=="3prong0pi0":
    #    varhelicity_plus_lab = ev.cosPsi_plus
    #if decaymodes[decaymode_p] in ['1prong1pi0', '1prong2pi0']:
    #    varhelicity_plus_lab = ev.energyasym_plus_lab

    #varhelicity_minus_lab = ev.evis_minus_lab
    #if decaymodes[decaymode_m]=="3prong0pi0":
    #    varhelicity_minus_lab = ev.cosPsi_minus
    #if decaymodes[decaymode_m] in ['1prong1pi0', '1prong2pi0']:
    #    varhelicity_minus_lab = ev.energyasym_minus_lab
    
    h_evis_corr[0][imode].Fill(ev.evis_minus_lab, ev.evis_plus_lab)
    h_energyasym_corr[0][imode].Fill(ev.energyasym_minus_lab, ev.energyasym_plus_lab)
    
    ## boson rest frame
    h_evis_corr[1][imode].Fill(ev.evis_minus_bRF, ev.evis_plus_bRF)
    h_energyasym_corr[1][imode].Fill(ev.energyasym_minus_bRF, ev.energyasym_plus_bRF)

    ## visible tau pair rest frame
    h_evis_corr[2][imode].Fill(ev.evis_minus_visRF, ev.evis_plus_visRF)
    h_energyasym_corr[2][imode].Fill(ev.energyasym_minus_visRF, ev.energyasym_plus_visRF)

    # costheta
    h_costheta_corr[imode].Fill(ev.cos_minus, ev.cos_plus)

    # cosPsi
    if combModes[imode] == "3prong0pi0_3prong0pi0":
        h_cosPsi_corr.Fill(ev.cosPsi_minus, ev.cosPsi_plus)

    # Old DecayMode
    index2_oDM = len(combModes)-1
    if decaymode_p in [1,2,4] and decaymode_m in [1,2,4]:
        ## lab frame
        h_evis_corr[0][index2_oDM].Fill(ev.evis_minus_lab, ev.evis_plus_lab)
        h_energyasym_corr[0][index2_oDM].Fill(ev.energyasym_minus_lab, ev.energyasym_plus_lab)
        ## boson rest frame
        h_evis_corr[1][index2_oDM].Fill(ev.evis_minus_bRF, ev.evis_plus_bRF)
        h_energyasym_corr[1][index2_oDM].Fill(ev.energyasym_minus_bRF, ev.energyasym_plus_bRF)
        ## visible tau pair rest frame
        h_evis_corr[2][index2_oDM].Fill(ev.evis_minus_visRF, ev.evis_plus_visRF)
        h_energyasym_corr[2][index2_oDM].Fill(ev.energyasym_minus_visRF, ev.energyasym_plus_visRF)

        h_costheta_corr[index2_oDM].Fill(ev.cos_minus, ev.cos_plus)

    ################
    # h_tauHelTag_prod[frame][mode]
    # h2d_tauHelTag_correlation[frame][mode]

    # lab frame
    h2d_tauHelTag_correlation[0][imode].Fill(heltagminus_lab, heltagplus_lab)
    h_tauHelTag_prod[0][imode].Fill(heltagminus_lab*heltagplus_lab)

    # boson rest frame
    h2d_tauHelTag_correlation[1][imode].Fill(heltagminus_bRF, heltagplus_bRF)
    h_tauHelTag_prod[1][imode].Fill(heltagminus_bRF*heltagplus_bRF)
    
    # Old DecayMode
    if decaymode_p in [1,2,4] and decaymode_m in [1,2,4]:
        # lab frame
        h2d_tauHelTag_correlation[0][index2_oDM].Fill(heltagminus_lab, heltagplus_lab)
        h_tauHelTag_prod[0][index2_oDM].Fill(heltagminus_lab*heltagplus_lab)
        # boson rest frame
        h2d_tauHelTag_correlation[1][index2_oDM].Fill(heltagminus_bRF, heltagplus_bRF)
        h_tauHelTag_prod[1][index2_oDM].Fill(heltagminus_bRF*heltagplus_bRF)
        
    #break

################
# draw and save histogmras

gStyle.SetOptStat(10)

canvas = TCanvas()

# 1D
for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for d, mode in enumerate(decaymodes):
            h_evis[f][c][d].Draw()
            name = h_evis[f][c][d].GetName().replace('h_','')
            canvas.SaveAs(args.outdir+name+'.png');
            h_energyasym[f][c][d].Draw()
            name = h_energyasym[f][c][d].GetName().replace('h_','')
            canvas.SaveAs(args.outdir+name+'.png')
for c, charge in enumerate(charges):
    for d, mode in enumerate(decaymodes):
        h_cos[c][d].Draw()
        name = h_cos[c][d].GetName().replace('h_','')
        canvas.SaveAs(args.outdir+name+'.png')
        
    h_cosPsi_3prongs[c].Draw()
    name = h_cosPsi_3prongs[c].GetName().replace('h_','')
    canvas.SaveAs(args.outdir+name+'.png')
        
# 2D
for m, mode in enumerate(combModes):
    h_costheta_corr[m].Draw('colz')
    name = h_costheta_corr[m].GetName().replace('h_','')
    canvas.SaveAs(args.outdir+name+'.png')

for f, frame in enumerate(frames):
    for m, mode in enumerate(combModes):
        h_evis_corr[f][m].Draw('colz')
        name = h_evis_corr[f][m].GetName().replace('h_','')
        canvas.SaveAs(args.outdir+name+'.png')
        h_energyasym_corr[f][m].Draw('colz')
        name = h_energyasym_corr[f][m].GetName().replace('h_','')
        canvas.SaveAs(args.outdir+name+'.png')
    
h_cosPsi_corr.Draw('colz')
name = h_cosPsi_corr.GetName().replace('h_','')
canvas.SaveAs(args.outdir+name+'.png')

for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        h_xThreeProngs[f][c].Draw('colz')
        name = h_xThreeProngs[f][c].GetName().replace('h_','')
        canvas.SaveAs(args.outdir+name+'.png')

##################
for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for d, mode in enumerate(decaymodes):
            h_tauHelTag[f][c][d].Draw()
            name = h_tauHelTag[f][c][d].GetName().replace('h_','')
            canvas.SaveAs(args.outdir+name+'.png')
        
for f, frame in enumerate(frames):
    for m, mode in enumerate(combModes):
        h2d_tauHelTag_correlation[f][m].Draw('colz')
        name = h2d_tauHelTag_correlation[f][m].GetName().replace('h2d_','')
        canvas.SaveAs(args.outdir+name+'.png')

        h_tauHelTag_prod[f][m].Draw('colz')
        name = h_tauHelTag_prod[f][m].GetName().replace('h_','')
        canvas.SaveAs(args.outdir+name+'.png')
