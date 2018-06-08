# FWLite script to print LHE header info from Marco Peruzzi

import ROOT
from DataFormats.FWLite import Events, Handle, Runs
import sys

#input_files = [
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/02259D47-A8D1-E611-AED1-02163E019BED.root', # 2016 sample
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/THQ_4f_Hincl_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/38553354-592A-E811-9526-A4BF011257F0.root', # 2017 sample
#]

input_files = sys.argv[1:]

is2016 = False

for myfile in input_files:
    events = Events(myfile)
    runs = Runs(myfile)
    lherunhandle = Handle('LHERunInfoProduct')
    for run in runs:
        run.getByLabel('externalLHEProducer' if not is2016 else 'source', lherunhandle)
        lheprod = lherunhandle.product()
        print '%d headers found'%lheprod.headers_size()
        it = lheprod.headers_begin()
        itidx=0
        while True:
            if it==lheprod.headers_end(): break
            print 'Printing tag #%d'%itidx,it.tag()
            for i in xrange(it.lines().size()):
                print it.lines().at(i)
            itidx+=1
            it+=1        
        print '%d comments found'%lheprod.comments_size()
        it = lheprod.comments_begin()
        itidx=0
        while True:
            if it==lheprod.comments_end(): break
            print 'Printing comment #%d'%itidx
            print it
        break
    genhandle = Handle('GenEventInfoProduct')
    lhehandle = Handle('LHEEventProduct')
    for event in events:
        event.getByLabel('generator', genhandle)
        event.getByLabel('externalLHEProducer' if not is2016 else 'source', lhehandle)
        genprod = genhandle.product()
        lheprod = lhehandle.product()
        print 'Weight in GenEventInfoProduct is',genprod.weight()
        ws = lheprod.weights()
        print 'Printing LHE weights: index/id/weight'
        for i in xrange(ws.size()):
            print i,ws.at(i).id, ws.at(i).wgt
        break
