#ifndef MVAEvaluator_h
#define MVAEvaluator_h

#include "TMVA/Reader.h"
#include "TString.h"

#include <string>

class MVAEvaluator
{
 public:
	MVAEvaluator();
	~MVAEvaluator(){};

	void setup_tmva_reader_HTT();
	void setup_tmva_reader_1l2tau_BDT1();
	void setup_tmva_reader_1l2tau_BDT2();
	void setup_tmva_reader_2lss1tau_BDT1();
	void setup_tmva_reader_2lss1tau_BDT2();
	void setup_tmva_reader_2lss1tau_BDT3();
	void setup_tmva_reader_2lss1tau_BDT4();
	void setup_tmva_reader_2lss1tau_BDT5();
	void setup_tmva_reader_2lss1tau_BDT6();
	void setup_tmva_reader_3l1tau_BDT1();
	void setup_tmva_reader_3l1tau_BDT2();
	void setup_tmva_reader_3l1tau_BDT3();
	void setup_tmva_reader_3l1tau_BDT4();
	void setup_tmva_reader_2l2tau_BDT1();
	void setup_tmva_reader_2l2tau_BDT2();
	void setup_tmva_reader_2l2tau_BDT3();
	void setup_tmva_reader_2l2tau_BDT4();

    float evaluate_bdt_HTT(float*);
	float evaluate_bdt_1l2tau_BDT1(float*);
	float evaluate_bdt_1l2tau_BDT2(float*);
	float evaluate_bdt_2lss1tau_BDT1(float*);
	float evaluate_bdt_2lss1tau_BDT2(float*);
	float evaluate_bdt_2lss1tau_BDT3(float*);
	float evaluate_bdt_2lss1tau_BDT4(float*);
	float evaluate_bdt_2lss1tau_BDT5(float*);
	float evaluate_bdt_2lss1tau_BDT6(float*);
	float evaluate_bdt_3l1tau_BDT1(float*);
	float evaluate_bdt_3l1tau_BDT2(float*);
	float evaluate_bdt_3l1tau_BDT3(float*);
	float evaluate_bdt_3l1tau_BDT4(float*);
	float evaluate_bdt_2l2tau_BDT1(float*);
	float evaluate_bdt_2l2tau_BDT2(float*);
	float evaluate_bdt_2l2tau_BDT3(float*);
	float evaluate_bdt_2l2tau_BDT4(float*);

 protected:

	static const TString data_directory_;
	
	float inputVars_HTT_[7];
	float inputVars_1l2tau_BDT1_[13];
	float inputVars_1l2tau_BDT2_[17];
	float inputVars_2lss1tau_BDT1_[15];
	float inputVars_2lss1tau_BDT2_[16];
	float inputVars_2lss1tau_BDT3_[18];
	float inputVars_2lss1tau_BDT4_[19];
	float inputVars_2lss1tau_BDT5_[20];
	float inputVars_2lss1tau_BDT6_[2];
	float inputVars_3l1tau_BDT1_[13];
	float inputVars_3l1tau_BDT2_[15];
	float inputVars_3l1tau_BDT3_[12];
	float inputVars_3l1tau_BDT4_[2];
	float inputVars_2l2tau_BDT1_[14];
	float inputVars_2l2tau_BDT2_[11];
	float inputVars_2l2tau_BDT3_[13];
	float inputVars_2l2tau_BDT4_[2];

	TMVA::Reader *reader_HTT_;
	TMVA::Reader *reader_1l2tau_BDT1_;
	TMVA::Reader *reader_1l2tau_BDT2_;
	TMVA::Reader *reader_2lss1tau_BDT1_;
	TMVA::Reader *reader_2lss1tau_BDT2_;
	TMVA::Reader *reader_2lss1tau_BDT3_;
	TMVA::Reader *reader_2lss1tau_BDT4_;
	TMVA::Reader *reader_2lss1tau_BDT5_;
	TMVA::Reader *reader_2lss1tau_BDT6_;
	TMVA::Reader *reader_3l1tau_BDT1_;
	TMVA::Reader *reader_3l1tau_BDT2_;
	TMVA::Reader *reader_3l1tau_BDT3_;
	TMVA::Reader *reader_3l1tau_BDT4_;
	TMVA::Reader *reader_2l2tau_BDT1_;
	TMVA::Reader *reader_2l2tau_BDT2_;
	TMVA::Reader *reader_2l2tau_BDT3_;
	TMVA::Reader *reader_2l2tau_BDT4_;
};

#endif
