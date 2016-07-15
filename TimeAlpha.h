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

// Enumerate variables for detector and model selection
enum tGetype { enrBEGe, enrCoax, natCoax, all };
enum tModel { LinAndExp, OnlyExp, OnlyLin };

// This is a TimeAlpha header file.
// Model source code is located in file TimeAlpha/TimeAlpha.cxx

// ---------------------------------------------------------
class TimeAlpha : public BCModel {
 public:

	// Constructor and destructor
	TimeAlpha( const char * name = "TimeAlpha" );
	~TimeAlpha();

	int ReadData( std::string keylist );

	void SetNBinsHistograms( int n, double min, double max )
	{
		fNBins = n;
		fHMinimum = min;
		fHMaximum = max;
	};

	void SetDetectorType( tGetype det ){ fDetectorType = det; };

	void SetModel( tModel model ){ fModel = model; };

	void SetNDetectors( int n ){ fNDetectors = n; };
	int GetNDetectors(){ return fNDetectors; };

	// Methods to overload, see file TimeAlpha.cxx
	double LogLikelihood(const std::vector<double> & parameters);

	// double LogAPrioriProbability(const std::vector<double> & parameters);

	void WriteOutput( std::string outputfilename );
	void WriteDistributions();

 private:

	tGetype fDetectorType;
	tModel fModel;

	int fNDetectors = 40;

	int fNBins;
	double fHMaximum, fHMinimum;

	// Data arrays with time distribution of events and livetime fraction in time (counting pulser events)
	TH1D * fHTimeAlpha;
	TH1D * fHLiveTimeFraction;

	std::vector<double> fVTimeAlpha;
	std::vector<double> fVLiveTimeFraction;

	int FillDataArrays();

};
// ---------------------------------------------------------

#endif
