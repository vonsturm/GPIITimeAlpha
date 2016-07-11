// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__TIMEALPHA__H
#define __BAT__TIMEALPHA__H

#include <cstdlib>
#include <vector>

#include "TH1D.h"

#include <BAT/BCModel.h>

// This is a TimeAlpha header file.
// Model source code is located in file TimeAlpha/TimeAlpha.cxx

// ---------------------------------------------------------
class TimeAlpha : public BCModel {
 public:

	// Constructor and destructor
	TimeAlpha( const char * name = "TimeAlpha", std::string keylist = "default.txt" );
	~TimeAlpha();

	int ReadData( std::string keylist );

	void SetNDetectors( int n ){ fNDetectors = n; };
	int GetNDetectors(){ return fNDetectors; };

	// Methods to overload, see file TimeAlpha.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

 private:

	int fNDetectors = 40;

	// Data arrays with time distribution of events and livetime fraction in time (counting pulser events)
	TH1D * fHTimeAlpha;
	TH1D * fHLiveTimeFraction;
	TH1D * fHRealDecay;

	std::vector<double> fVTimeAlpha;
	std::vector<double> fVLiveTimeFraction;
	std::vector<double> fVRealDecay;

};
// ---------------------------------------------------------

#endif
