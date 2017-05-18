import numpy as np
#from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from ttHTauTauAnalysis.ttHtautauAnalyzer.mva_utils import get_inputs, plot_roc, plot_clf_results, plot_correlation, print_variable_rank

variables = """
nJet
mindr_lep0_jet
mindr_lep1_jet
avg_dr_jet
max_lep_eta
met
mht
mT_met_lep0
lep0_conept
lep1_conept
dr_leps
dr_lep0_tau
dr_lep1_tau
""".split()

xsig, ysig, wsig = get_inputs('ttH', variables)
xbkg, ybkg, wbkg = get_inputs('ttV', variables)

plot_correlation(xsig, variables, 'correlations_sig.png')
plot_correlation(xbkg, variables, 'correlations_bkg.png')

x = np.concatenate((xsig, xbkg))
y = np.concatenate((ysig, ybkg))
w = np.concatenate((wsig, wbkg))

x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(x, y, w, test_size=0.2)

bdt = GradientBoostingClassifier(n_estimators=100, learning_rate=0.2, max_depth=3, random_state=0)
bdt.fit(x_train, y_train)#, w_train)

joblib.dump(bdt, 'bdt.pkl')
#bdt=joblib.load('bdt.pkl')

# evaluate training results
#y_decision = bdt.decision_function(x_test)#.ravel()
y_proba_sig = bdt.predict_proba(x_test)[:,0]
plot_clf_results(bdt, x_train, y_train, w_train, x_test, y_test, w_test)
#plot_roc(y_test, y_decision, w_test, 'roc.png')
plot_roc(y_test, y_proba_sig, w_test, 'roc.png')
print_variable_rank(bdt, variables)
