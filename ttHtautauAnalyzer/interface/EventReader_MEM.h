#ifndef EventReader_MEM_h
#define EventReader_MEM_h

//#include "ttH_Htautau_MEM_Analysis/MEM/interface/PyRun2EventData_t.h"
//#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/NodeScheduler.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/MGIntegration.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/eventNtuple.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <string>
#include <algorithm>

class EventReader_MEM
{
 public:
	
	EventReader_MEM(const TString&, int64_t, int64_t, bool, const TString&);
	~EventReader_MEM();

	void fillEventList(IntegrationMsg_t*, int64_t);

 private:

	TFile* infile_;
	TTree* intree_;
	eventNtuple evNtuple_;
	
	int64_t startevent_;
	int64_t nevents_;

	bool looseSel_;
};

#endif
