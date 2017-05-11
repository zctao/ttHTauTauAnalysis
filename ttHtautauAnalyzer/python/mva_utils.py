import numpy as np
import pandas as pd
import ROOT as r
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from root_numpy import fill_hist
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import classification_report

def plot_correlation(data,figname,**kwds):
    """Calculate pairwise correlation between variables"""

    corrmat = data.corr(**kwds)

    fig, ax1 = plt.subplots(ncols=1, figsize=(6,5))

    opts = {'cmap': plt.get_cmap("RdBu"),'vmin': -1, 'vmax': +1}
    heatmap1 = ax1.pcolor(corrmat, **opts)
    plt.colorbar(heatmap1, ax=ax1)

    ax1.set_title("Correlations")

    labels = corrmat.columns.values
    for ax in (ax1,) :
        # shift location of ticks to center of the bins
        ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
        ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
        ax.set_xticklabels(labels, minor=False, ha='right', rotation=70)
        ax.set_yticklabels(labels, minor=False)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def plot_roc(y_test, y_decision, figname):
    fpr, tpr, thresholds = roc_curve(y_test, y_decision)
    roc_auc = roc_auc_score(y_test, y_decision)
    
    plt.plot(fpr, tpr, lw=1, label='ROC (area = %.02f)'%(roc_auc))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Background efficiency')
    plt.ylabel('Signal efficiency')
    plt.legend(loc="lower right")
    plt.grid()
    #plt.show()
    plt.savefig(figname)
    plt.close()

def plot_clf_results(clf, x_train, y_train, w_train, x_test, y_test, w_test, nbins=30):
    decisions = []
    weights = []
    for x,y,w in ((x_train, y_train, w_train), (x_test, y_test, w_test)):
        dsig = clf.decision_function(x[y>0])
        wsig = w[y>0]
        dbkg = clf.decision_function(x[y<0])
        wbkg = w[y<0]
        decisions += [dsig,dbkg]
        weights += [wsig, wbkg]

    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    
    h_train_sig = r.TH1D('train_sig','',nbins,low,high)
    h_train_sig.SetStats(0)
    fill_hist(h_train_sig, decisions[0], weights=weights[0])

    h_train_bkg = r.TH1D('train_bkg','',nbins,low,high)
    h_train_bkg.SetStats(0)
    fill_hist(h_train_bkg, decisions[1], weights=weights[1])

    h_test_sig = r.TH1D('test_sig','',nbins,low,high)
    h_test_sig.SetStats(0)
    fill_hist(h_test_sig, decisions[2], weights=weights[2])

    h_test_bkg = r.TH1D('test_bkg','',nbins,low,high)
    h_test_bkg.SetStats(0)
    fill_hist(h_test_bkg, decisions[3], weights=weights[3])

    # legend
    l = r.TLegend(0.75,0.75,0.9,0.9)
    l.AddEntry(h_train_sig,'S (train)','f')
    l.AddEntry(h_train_bkg,'B (train)','f')
    l.AddEntry(h_test_sig,'S (test)','p')
    l.AddEntry(h_test_bkg,'B (test)','p')
    
    # draw histograms
    tcanvas = r.TCanvas()
    h_train_sig.SetLineColor(2)
    h_train_sig.SetFillColorAlpha(2,0.5)
    h_train_sig.GetXaxis().SetTitle('BDT Output')
    h_train_sig.GetYaxis().SetTitle('A.U.')
    h_train_sig.Draw("HIST")
    h_train_bkg.SetLineColor(4)
    h_train_bkg.SetFillColorAlpha(4,0.5)
    h_train_bkg.Draw('SAME HIST')
    h_test_sig.SetMarkerStyle(20)
    h_test_sig.SetMarkerColor(2)
    h_test_sig.SetMarkerSize(0.8)
    h_test_sig.SetLineColor(2)
    h_test_sig.Draw('SAME')
    h_test_bkg.SetMarkerStyle(20)
    h_test_bkg.SetMarkerColor(4)
    h_test_bkg.SetMarkerSize(0.8)
    h_test_bkg.SetLineColor(4)
    h_test_bkg.Draw('SAME')
    l.Draw('SAME')

    tcanvas.SaveAs('BDTOutput.png')
    
def evaluate_inputs(x, y, variables):
    
    df = pd.DataFrame(np.hstack((x, y.reshape(y.shape[0], -1))),
                  columns=variables+['y'])
    
    sig = df.y > 0.
    bkg = df.y < 0.
    
    plot_correlation(df[sig].drop('y',1),'correlations_sig.png')
    plot_correlation(df[bkg].drop('y',1),'correlations_bkg.png')
    
