import argparse
from ROOT import TCanvas, TFile, TTree, TH1F, TH2F
from ROOT import gStyle

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

# h_upsilon[frame][tau charge][decay mode]
h_upsilon = [[[] for j in range(len(charges))] for i in range(len(frames))]

for f, frame in enumerate(frames):
    for c, charge in enumerate(charges):
        for decay in decaymodes:
            h_tmp = TH1F('h_upsilon_'+charge+'_'+decay+"_"+frame,
                         '#tau^{'+signs[charge]+'}  '+decay+'  '+frame+' frame',
                         20, -1., 1.+args.epsilon)
            h_tmp.GetXaxis().SetTitle("(E_{#pm}-E_{0})/(E_{#pm}+E_{0})")
            h_upsilon[f][c].append(h_tmp)
            
#h_cos[tau charge][decay mode]
h_cos = [[] for i in range(len(charges))]

for c, charge in enumerate(charges):
    for decay in decaymodes:
        h_tmp = TH1F('h_cos_'+charge+'_'+decay,'#tau^{'+signs[charge]+'}  '+decay,
                     20, -1., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle("cos#theta^{"+signs[charge]+"}")
        h_cos[c].append(h_tmp)

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

# h_upsilon_corr[frame][mode]
h_upsilon_corr = [[] for i in range(len(frames))]
for f, frame in enumerate(frames):
    for mode in combModes:
        h_tmp = TH2F('h_upsilon_corr_'+mode+'_'+frame, mode+' ('+frame+' frame)',
                     20, -1., 1.+args.epsilon, 20, -1., 1.+args.epsilon)
        h_tmp.GetXaxis().SetTitle('#Upsilon^{-} ('+mode.split('_')[1]+')')
        h_tmp.GetYaxis().SetTitle('#Upsilon^{+} ('+mode.split('_')[0]+')')
        h_upsilon_corr[f].append(h_tmp)

# h_evis_upsilon_corr[frame][mode]
h_evis_upsilon_corr = [[] for i in range(len(frames))]
for f, frame in enumerate(frames):
    for mode in combModes:
        xmin = 0.
        xmax = 50. if frame=='lab' else 1.+args.epsilon
        ymin = 0.
        ymax = 50. if frame=='lab' else 1.+args.epsilon
        xtitle = '2E_{vis}^{-}/M_{X}'
        ytitle = '2E_{vis}^{+}/M_{X}'
        
        if mode.split('_')[1] in ['1prong1pi0','1prong2pi0']:
            xmin = -1.
            xmax = 1.+args.epsilon
            xtitle = '#Upsilon^{-}'
        if mode.split('_')[0] in ['1prong1pi0','1prong2pi0']:
            ymin = -1.
            ymax = 1.+args.epsilon
            ytitle = '#Upsilon^{+}'
            
        h_tmp = TH2F('h_evis_upsilon_corr_'+mode+'_'+frame,
                     mode+' ('+frame+' frame)', 20, xmin, xmax, 20, ymin, ymax)
        h_tmp.GetXaxis().SetTitle(xtitle+' ('+mode.split('_')[1]+')')
        h_tmp.GetYaxis().SetTitle(ytitle+' ('+mode.split('_')[0]+')')
        h_evis_upsilon_corr[f].append(h_tmp)
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
    h_upsilon[0][0][decaymode_p].Fill(ev.upsilon_plus_lab)
    # tau-
    h_evis[0][1][decaymode_m].Fill(ev.evis_minus_lab)
    h_upsilon[0][1][decaymode_m].Fill(ev.upsilon_minus_lab)
    
    ## boson rest frame
    # tau+
    h_evis[1][0][decaymode_p].Fill(ev.evis_plus_bRF)
    h_upsilon[1][0][decaymode_p].Fill(ev.upsilon_plus_bRF)
    # tau-
    h_evis[1][1][decaymode_m].Fill(ev.evis_minus_bRF)
    h_upsilon[1][1][decaymode_m].Fill(ev.upsilon_minus_bRF)
    
    ## visible tau pair rest frame
    # tau+
    h_evis[2][0][decaymode_p].Fill(ev.evis_plus_visRF)
    h_upsilon[2][0][decaymode_p].Fill(ev.upsilon_plus_visRF)
    # tau-
    h_evis[2][1][decaymode_m].Fill(ev.evis_minus_visRF)
    h_upsilon[2][1][decaymode_m].Fill(ev.upsilon_minus_visRF)
    
    h_cos[0][decaymode_p].Fill(ev.cos_plus)
    h_cos[1][decaymode_m].Fill(ev.cos_minus)

    # Old DecayMode
    index_oldDM = len(decaymodes)-1
    if decaymode_p in [1,2,4] and decaymode_m in [1,2,4]:
        h_evis[0][0][index_oldDM].Fill(ev.evis_plus_lab)
        h_upsilon[0][0][index_oldDM].Fill(ev.upsilon_plus_lab)
        h_evis[1][0][index_oldDM].Fill(ev.evis_plus_bRF)
        h_upsilon[1][0][index_oldDM].Fill(ev.upsilon_plus_bRF)
        h_evis[2][0][index_oldDM].Fill(ev.evis_plus_visRF)
        h_upsilon[2][0][index_oldDM].Fill(ev.upsilon_plus_visRF)
        
        h_evis[0][1][index_oldDM].Fill(ev.evis_plus_lab)
        h_upsilon[0][1][index_oldDM].Fill(ev.upsilon_plus_lab)
        h_evis[1][1][index_oldDM].Fill(ev.evis_minus_bRF)
        h_upsilon[1][1][index_oldDM].Fill(ev.upsilon_minus_bRF)
        h_evis[2][1][index_oldDM].Fill(ev.evis_minus_visRF)
        h_upsilon[2][1][index_oldDM].Fill(ev.upsilon_minus_visRF)

    ### Correlations    
    if decaymode_p > decaymode_m:
        continue  # only plot e.g. l+pi-, but not l-pi+
    
    # Index of tau pair decay mode in combModes list
    #(2*len(combModes)-decaymode_p+1)*decaymode_p/2 + decaymode_m - decaymode_p
    imode = (2*(len(decaymodes)-1)-decaymode_p-1) * decaymode_p / 2 + decaymode_m
    #print decaymodes[decaymode_p], decaymodes[decaymode_m]
    #print combModes[imode]
    assert(decaymodes[decaymode_p]+'_'+decaymodes[decaymode_m]==combModes[imode])
    plotUpsilon_x = decaymodes[decaymode_m] in ['1prong1pi0', '1prong2pi0']
    assert(not(plotUpsilon_x ^ ('Upsilon' in h_evis_upsilon_corr[0][imode].GetXaxis().GetTitle())))
    plotUpsilon_y = decaymodes[decaymode_p] in ['1prong1pi0', '1prong2pi0']
    assert(not(plotUpsilon_y ^ ('Upsilon' in h_evis_upsilon_corr[0][imode].GetYaxis().GetTitle())))
    
    ## lab frame
    h_evis_corr[0][imode].Fill(ev.evis_minus_lab, ev.evis_plus_lab)
    h_upsilon_corr[0][imode].Fill(ev.upsilon_minus_lab, ev.upsilon_plus_lab)
    h_evis_upsilon_corr[0][imode].Fill(
        ev.upsilon_minus_lab if plotUpsilon_x else ev.evis_minus_lab,
        ev.upsilon_plus_lab if plotUpsilon_y else ev.evis_plus_lab)
    
    ## boson rest frame
    h_evis_corr[1][imode].Fill(ev.evis_minus_bRF, ev.evis_plus_bRF)
    h_upsilon_corr[1][imode].Fill(ev.upsilon_minus_bRF, ev.upsilon_plus_bRF)
    h_evis_upsilon_corr[1][imode].Fill(
        ev.upsilon_minus_bRF if plotUpsilon_x else ev.evis_minus_bRF,
        ev.upsilon_plus_bRF if plotUpsilon_y else ev.evis_plus_bRF)

    ## visible tau pair rest frame
    h_evis_corr[2][imode].Fill(ev.evis_minus_visRF, ev.evis_plus_visRF)
    h_upsilon_corr[2][imode].Fill(ev.upsilon_minus_visRF, ev.upsilon_plus_visRF)
    h_evis_upsilon_corr[2][imode].Fill(
        ev.upsilon_minus_visRF if plotUpsilon_x else ev.evis_minus_visRF,
        ev.upsilon_plus_visRF if plotUpsilon_y else ev.evis_plus_visRF)


    # costheta
    h_costheta_corr[imode].Fill(ev.cos_minus, ev.cos_plus)

    # Old DecayMode
    index2_oDM = len(combModes)-1
    if decaymode_p in [1,2,4] and decaymode_m in [1,2,4]:
        ## lab frame
        h_evis_corr[0][index2_oDM].Fill(ev.evis_minus_lab, ev.evis_plus_lab)
        h_upsilon_corr[0][index2_oDM].Fill(ev.upsilon_minus_lab, ev.upsilon_plus_lab)
        h_evis_upsilon_corr[0][index2_oDM].Fill(
            ev.upsilon_minus_lab if plotUpsilon_x else ev.evis_minus_lab,
            ev.upsilon_plus_lab if plotUpsilon_y else ev.evis_plus_lab)
        ## boson rest frame
        h_evis_corr[1][index2_oDM].Fill(ev.evis_minus_bRF, ev.evis_plus_bRF)
        h_upsilon_corr[1][index2_oDM].Fill(ev.upsilon_minus_bRF, ev.upsilon_plus_bRF)
        h_evis_upsilon_corr[1][index2_oDM].Fill(
            ev.upsilon_minus_bRF if plotUpsilon_x else ev.evis_minus_bRF,
            ev.upsilon_plus_bRF if plotUpsilon_y else ev.evis_plus_bRF)
        ## visible tau pair rest frame
        h_evis_corr[2][index2_oDM].Fill(ev.evis_minus_visRF, ev.evis_plus_visRF)
        h_upsilon_corr[2][index2_oDM].Fill(ev.upsilon_minus_visRF, ev.upsilon_plus_visRF)
        h_evis_upsilon_corr[2][index2_oDM].Fill(
            ev.upsilon_minus_visRF if plotUpsilon_x else ev.evis_minus_visRF,
            ev.upsilon_plus_visRF if plotUpsilon_y else ev.evis_plus_visRF)

        h_costheta_corr[index2_oDM].Fill(ev.cos_minus, ev.cos_plus)
    
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
            name = h_evis[f][c][d].GetName().lstrip('h_')
            canvas.SaveAs(args.outdir+name+'.png');
            h_upsilon[f][c][d].Draw()
            name = h_upsilon[f][c][d].GetName().lstrip('h_')
            canvas.SaveAs(args.outdir+name+'.png')
for c, charge in enumerate(charges):
    for d, mode in enumerate(decaymodes):
        h_cos[c][d].Draw()
        name = h_cos[c][d].GetName().lstrip('h_')
        canvas.SaveAs(args.outdir+name+'.png')

# 2D
for m, mode in enumerate(combModes):
    h_costheta_corr[m].Draw('colz')
    name = h_costheta_corr[m].GetName().lstrip('h_')
    canvas.SaveAs(args.outdir+name+'.png')

for f, frame in enumerate(frames):
    for m, mode in enumerate(combModes):
        h_evis_corr[f][m].Draw('colz')
        name = h_evis_corr[f][m].GetName().lstrip('h_')
        canvas.SaveAs(args.outdir+name+'.png')
        h_upsilon_corr[f][m].Draw('colz')
        name = h_upsilon_corr[f][m].GetName().lstrip('h_')
        canvas.SaveAs(args.outdir+name+'.png')
        h_evis_upsilon_corr[f][m].Draw('colz')
        name = h_evis_upsilon_corr[f][m].GetName().lstrip('h_')
        canvas.SaveAs(args.outdir+name+'.png')
    
