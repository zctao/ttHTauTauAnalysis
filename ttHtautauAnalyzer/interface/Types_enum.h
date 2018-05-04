#ifndef Types_enum_h
#define Types_enum_h

/*
 * enum for analysis type.
 */

enum Analysis_types {
	Analyze_lepton_jet,
	Analyze_dilepton,
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
	Signal_2lss1tau,
	Signal_1l2tau,
	Signal_3l1tau,
	Signal_2l2tau,
	Control_2los1tau,
	Control_fake_2lss1tau,
	Control_fake_1l2tau,
	Control_fake_3l1tau,
	Control_fake_2l2tau,
	Control_WZ,
	Loose_2lss1tau,  // for training
	Loose_1l2tau,
	Loose_3l1tau,
	Inclusive_1l2tau,
	Inclusive_2lss1tau,
	Inclusive_3l1tau,
	Inclusive_2l2tau,
    Control_ttW,
	Control_ttZ,
	Selection_NA
};


#endif
