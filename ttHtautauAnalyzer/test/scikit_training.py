import numpy as np
from root_numpy import root2array, rec2array
#from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from ttHTauTauAnalysis.ttHtautauAnalyzer.mva_utils import plot_roc, plot_clf_results , evaluate_inputs

fn_sig = "mvaVars_ttH_loose.root"
fn_bkg = "mvaVars_TTZ_loose.root"

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

array_sig = rec2array(root2array(fn_sig, 'mva', variables))
array_wsig = root2array(fn_sig, 'mva', 'event_weight')
array_bkg = rec2array(root2array(fn_bkg, 'mva', variables))
array_wbkg = root2array(fn_bkg, 'mva', 'event_weight')

# normalize weights
array_wsig *= 1./np.sum(array_wsig)
array_wbkg *= 1./np.sum(array_wbkg)

x = np.concatenate((array_sig, array_bkg))
y = np.concatenate([np.ones(array_sig.shape[0]), -1*np.ones(array_bkg.shape[0])])
w = np.concatenate((array_wsig, array_wbkg))

evaluate_inputs(x, y, variables)

x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(x, y, w, test_size=0.5)

bdt = GradientBoostingClassifier(n_estimators=100, learning_rate=0.2, max_depth=3, random_state=0)

bdt.fit(x_train, y_train)#, w_train)

joblib.dump(bdt, 'bdt.pkl')

# evaluate training results
y_decision = bdt.decision_function(x_test)#.ravel()
plot_clf_results(bdt, x_train, y_train, w_train, x_test, y_test, w_test)
plot_roc(y_test, y_decision,'roc.png')
