import numpy as np
import pandas as pd
import ROOT as r
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from root_numpy import root2array, rec2array, fill_hist
from sklearn.metrics import roc_curve, roc_auc_score, auc
from sklearn.metrics import classification_report
from sklearn.model_selection import learning_curve
from sklearn.model_selection import GridSearchCV

# mostly based on:
# https://betatim.github.io/posts/sklearn-for-TMVA-users/
# https://betatim.github.io/posts/advanced-sklearn-for-TMVA/

def get_inputs(sample_name, variables, tree_name='mva'):
    x = None
    y = None
    w = None

    infiles = []
    xsections = []
    if ('ttH' in sample_name):
        infiles = ["mvaVars_ttH_loose.root"]
        xsections = [0.215]
    elif ('ttV' in sample_name):
        infiles = ["mvaVars_TTZ_loose.root", "mvaVars_TTW_loose.root"]
        xsections = [0.253, 0.204]  # [TTZ, TTW]
    elif ('ttbar' in sample_name):
        infiles = ["mvaVars_TTSemilep_loose.root","mvaVars_TTDilep_loose.root"]
        xsections = [182, 87.3]  # [semilep, dilep]
    else:
        print "Pick one sample name from 'ttH', 'ttV' or 'ttbar'"
        return x, y, w

    for fn, xs in zip(infiles, xsections):
        xi = rec2array(root2array(fn, tree_name, variables))
        wi = root2array(fn, tree_name, 'event_weight')

        # scale weight and renormalize total weights to one
        wi *= (xs / sum(xsections)) /np.sum(wi)

        if x is not None:
            x = np.concatenate((x,xi))
            w = np.concatenate((w,wi))
        else:
            x = xi
            w = wi

    if ('ttH' in sample_name):
        y = np.ones(x.shape[0])
    else:
        y = np.zeros(x.shape[0])
        #y = -1*np.ones(x.shape[0])

    return x, y, w


def plot_correlation(x, variables, figname, **kwds):
    df = pd.DataFrame(x, columns=variables)
    
    """Calculate pairwise correlation between variables"""
    corrmat = df.corr(**kwds)

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
    
    
def plot_roc(y_test, y_pred, w_test, figname):
    fpr, tpr, thresholds = roc_curve(y_test, y_pred, sample_weight=w_test)
    print max(thresholds)
    print min(thresholds)
    roc_auc = roc_auc_score(y_test, y_pred, sample_weight=w_test)
    #roc_auc = auc(fpr, tpr,reorder=True)
    
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
        w *= 1./np.sum(w)
        #dsig = clf.decision_function(x[y>0.5])
        dsig = clf.predict_proba(x[y>0.5])[:,0]
        wsig = w[y>0.5]
        #dbkg = clf.decision_function(x[y<0.5])
        dbkg = clf.predict_proba(x[y<0.5])[:,0]
        wbkg = w[y<0.5]
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

    
def print_variable_rank(clf, variables):
    print 'classifier variable ranking : '
    for var, score in sorted(zip(variables, clf.feature_importances_),key=lambda x: x[1], reverse=True):
        print var, score

        
def run_grid_search(clf, x_dev, y_dev, verbose=True):
    param_grid = {"n_estimators": [50,200,400,1000],
              "max_depth": [1, 3, 8],
              'learning_rate': [0.1, 0.2, 1.]}
    clfGS = GridSearchCV(clf, param_grid, cv=3, scoring='roc_auc',
                                    n_jobs=8)
    clfGS.fit(x_dev,y_dev)

    print 'Best set of parameters : '
    print clfGS.best_estimator_

    if verbose:
        print
        print 'Grid scores on a subset of the development set:'
        for params, mean_score, scores in clfGS.grid_scores_:
            print "%0.4f (+/-%0.04f) for %r"%(mean_score, scores.std(), params)
        
        #y_true, y_pred = y_dev, clf.decision_function(x_dev)
        #print "  It scores %0.4f on the full development set"%roc_auc_score(y_true, y_pred)
        #y_true, y_pred = y_eval, clf.decision_function(x_eval)
        #print "  It scores %0.4f on the full evaluation set"%roc_auc_score(y_true, y_pred)

        
def plot_validation_curve(clfs, x_train, y_train, x_test, y_test,
                          figname="validation_curve.png"):
    for n,clf in enumerate(clfs):
        test_score = np.empty(len(clf.estimators_))
        train_score = np.empty(len(clf.estimators_))
        
        for i, pred in enumerate(clf.staged_decision_function(x_test)):
            test_score[i] = roc_auc_score(y_test, pred)

        for i, pred in enumerate(clf.staged_decision_function(x_train)):
            train_score[i] = roc_auc_score(y_train, pred)

        best_iter = np.argmax(test_score)
        rate = clf.get_params()['learning_rate']
        depth = clf.get_params()['max_depth']
        test_line = plt.plot(test_score,label='rate=%.1f depth=%i (%.2f)'%(rate,depth,test_score[best_iter]))
        colour = test_line[-1].get_color()
        plt.plot(train_score, '--', color=colour)
        plt.xlabel("Number of boosting iterations")
        plt.ylabel("Area under ROC")
        plt.axvline(x=best_iter, color=colour)

    plt.legend(loc='best')
    return plt
    #plt.savefig(figname)
    #plt.close()

    
def plot_learning_curve(estimator, title, X, y, cv=None,
                        n_jobs=1, train_sizes=np.linspace(.1, 1.0, 5),
                        scoring=None, ax=None, xlabel=True):
    if ax is None:
        plt.figure()
        ax.title(title)
    
    if xlabel:
        ax.set_xlabel("Training examples")
        
    ax.set_ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(estimator,
                                                            X, y,
                                                            cv=cv,
                                                            n_jobs=n_jobs,
                                                            train_sizes=train_sizes,
                                                            scoring=scoring)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)

    ax.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    ax.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    ax.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    ax.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    ax.set_ylim([0.65, 1.0])
    return plt
