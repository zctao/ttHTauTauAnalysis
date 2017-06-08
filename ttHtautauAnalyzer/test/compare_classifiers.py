import numpy as np
import ttHTauTauAnalysis.ttHtautauAnalyzer.mva_utils as util

# sklearn weighted
dir = '/uscms/home/ztao/nobackup/BDTs/sklearn/ttV/weighted/'
name = 'sklearn_weighted'
result_sklearn_weighted = util.get_sklearn_test_results(dir, name)
# get event weights
y_test, y_pred, event_weights, label = result_sklearn_weighted

# sklearn unweighted
dir = '/uscms/home/ztao/nobackup/BDTs/sklearn/ttV/unweighted/'
name = 'sklearn_unweighted'
result_sklearn_unweighted = util.get_sklearn_test_results(dir, name)

# sklearn  flip negative weights
dir = '/uscms/home/ztao/nobackup/BDTs/sklearn/ttV/flipnegweight/'
name = 'sklearn_flip'
result_sklearn_flip = util.get_sklearn_test_results(dir, name)

# sklearn ignore negative weights
dir = '/uscms/home/ztao/nobackup/BDTs/sklearn/ttV/ignorenegweight/'
name = 'sklearn_ignore'
result_sklearn_ignore = util.get_sklearn_test_results(dir, name)

# tmva inverse boost negative weights
dir = '/uscms/home/ztao/nobackup/BDTs/tmva/ttV/weighted/'
name = 'tmva_weighted'
result_tmva_weighted = util.get_tmva_test_results(dir, util.variables_ttV, name)

# tmva unweighted
dir = '/uscms/home/ztao/nobackup/BDTs/tmva/ttV/unweighted/'
name = 'tmva_unweighted'
result_tmva_unweighted = util.get_tmva_test_results(dir, util.variables_ttV, name)

# flip negative weights
dir = '/uscms/home/ztao/nobackup/BDTs/tmva/ttV/flipnegweight/'
name = 'tmva_flip'
result_tmva_flip = util.get_tmva_test_results(dir, util.variables_ttV, name)

# ignore negative weights
dir = '/uscms/home/ztao/nobackup/BDTs/tmva/ttV/ignorenegweight/'
name = 'tmva_ignore'
result_tmva_ignore = util.get_tmva_test_results(dir, util.variables_ttV, name)

# pair annihilation
#dir = '/uscms/home/ztao/nobackup/BDTs/tmva/ttV/pairnegweight/'
#name = 'tmva_pair'
#result_tmva_pair = util.get_tmva_test_results(dir, util.variables_ttV, name)

#util.plot_rocs([result_sklearn_weighted, result_sklearn_unweighted, result_sklearn_flip, result_sklearn_ignore], 'weights_sklearn_roc.png')
#util.plot_rocs([result_tmva_weighted, result_tmva_unweighted, result_tmva_flip, result_tmva_ignore, result_tmva_pair], 'weights_tmva_roc.png')
util.plot_rocs([result_sklearn_weighted, result_sklearn_unweighted, result_sklearn_flip, result_sklearn_ignore, result_tmva_weighted, result_tmva_unweighted, result_tmva_flip, result_tmva_ignore], 'compare_rocs.png', weights=event_weights)
