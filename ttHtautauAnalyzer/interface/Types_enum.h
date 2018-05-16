#ifndef Types_enum_h
#define Types_enum_h

#include <map>
#include <string>

/*
 * enum for analysis type.
 */

enum Analysis_types {
	Analyze_1l2tau,
	Analyze_2lss1tau,
	Analyze_3l1tau,
	Analyze_2l2tau,
	Analyze_inclusive,
	Analyze_NA
};

/*
 * enum for event selection type.
 */

// enum for event selection region
enum Selection_types {
	Signal_1l2tau,
	Signal_2lss1tau,
	Signal_3l1tau,
	Signal_2l2tau,
	Application_Fake_1l2tau,
	Application_Fake_2lss1tau,
	Application_Flip_2lss1tau,
	Application_Fake_3l1tau,
	Application_Fake_2l2tau,
	Control_Fake_1l2tau,
	Control_Fake_2lss1tau,
	Control_Fake_3l1tau,
	Control_Fake_2l2tau,
	Control_FakeAR_1l2tau,
	Control_FakeAR_2lss1tau,
	Control_FakeAR_3l1tau,
	Control_FakeAR_2l2tau,
	Loose_2lss1tau,  // for training
	Loose_1l2tau,
	Loose_3l1tau,
	Loose_2l2tau,
	Inclusive_1l2tau,
	Inclusive_2lss1tau,
	Inclusive_3l1tau,
	Inclusive_2l2tau,
    Control_ttW,
	Control_ttZ,
	Control_WZ,
	Selection_NA
};

const std::map<std::string, Analysis_types> AnaTypeMap = {
	{"1l2tau", Analysis_types::Analyze_1l2tau},
	{"2lss1tau", Analysis_types::Analyze_2lss1tau},
	{"3l1tau", Analysis_types::Analyze_3l1tau},
	{"2l2tau", Analysis_types::Analyze_2l2tau},
	{"inclusive", Analysis_types::Analyze_inclusive},
	{"NA", Analysis_types::Analyze_NA}
};

const std::map<std::string, Selection_types> SelTypeMap = {
	{"signal_1l2tau", Selection_types::Signal_1l2tau},
	{"application_fake_1l2tau", Selection_types::Application_Fake_1l2tau},
	{"control_fake_1l2tau", Selection_types::Control_Fake_1l2tau},
	{"control_fakeAR_1l2tau", Selection_types::Control_FakeAR_1l2tau},
	{"signal_2lss1tau", Selection_types::Signal_2lss1tau},
	{"application_fake_2lss1tau", Selection_types::Application_Fake_2lss1tau},
	{"application_flip_2lss1tau", Selection_types::Application_Flip_2lss1tau},
	{"control_fake_2lss1tau", Selection_types::Control_Fake_2lss1tau},
	{"control_fakeAR_2lss1tau", Selection_types::Control_FakeAR_2lss1tau},
	{"signal_3l1tau", Selection_types::Signal_3l1tau},
	{"application_fake_3l1tau", Selection_types::Application_Fake_3l1tau},
	{"control_fake_3l1tau", Selection_types::Control_Fake_3l1tau},
	{"control_fakeAR_3l1tau", Selection_types::Control_FakeAR_3l1tau},
	{"signal_2l2tau", Selection_types::Signal_2l2tau},
	{"application_fake_2l2tau", Selection_types::Application_Fake_2l2tau},
	{"control_fake_2l2tau", Selection_types::Control_Fake_2l2tau},
	{"control_fakeAR_2l2tau", Selection_types::Control_FakeAR_2l2tau},
	{"loose_1l2tau", Selection_types::Loose_1l2tau},
	{"loose_2lss1tau", Selection_types::Loose_2lss1tau},
	{"loose_3l1tau", Selection_types::Loose_3l1tau},
	{"loose_2l2tau", Selection_types::Loose_2l2tau},
	{"inclusive_1l2tau", Selection_types::Inclusive_1l2tau},
	{"inclusive_2lss1tau", Selection_types::Inclusive_2lss1tau},
	{"inclusive_3l1tau", Selection_types::Inclusive_3l1tau},
	{"inclusive_2l2tau", Selection_types::Inclusive_2l2tau},
	{"control_ttW", Selection_types::Control_ttW},
	{"control_ttZ", Selection_types::Control_ttZ},
	{"control_WZ", Selection_types::Control_WZ},
	{"NA", Selection_types::Selection_NA}
};

namespace Types_enum {
	inline Analysis_types getAnaType(const std::string& anatype) {
		Analysis_types anaType = Analyze_NA;
		if (AnaTypeMap.count(anatype)>0) anaType = AnaTypeMap.at(anatype);
		return anaType;
	}

	inline Selection_types getSelType(const std::string& seltype) {
		Selection_types selType = Selection_NA;
		if (SelTypeMap.count(seltype)>0) selType = SelTypeMap.at(seltype);
		return selType;
	}

	inline int getNnominalLeptons(Analysis_types anatype) {
		int nleps = -1;
		if (anatype==Analyze_1l2tau)
			nleps = 1;
		else if (anatype==Analyze_2lss1tau or anatype==Analyze_2l2tau)
			nleps = 2;
		else if (anatype==Analyze_3l1tau)
			nleps = 3;

		return nleps;
	}

	inline int getNnominalTaus(Analysis_types anatype, Selection_types seltype) {
		int ntaus = -1;
		if (anatype==Analyze_1l2tau or anatype==Analyze_2l2tau)
			ntaus = 2;
		else if (anatype==Analyze_2lss1tau or anatype==Analyze_3l1tau) {
			if (seltype==Control_ttW or seltype==Control_ttZ or seltype==Control_WZ)
				ntaus = 0;
			else
				ntaus = 1;
		}
		return ntaus;
	}
}

#endif
