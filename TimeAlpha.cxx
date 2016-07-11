// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <string>

#include "TimeAlpha.h"

#include <BAT/BCMath.h>

#include "TChain.h"
#include "TH1D.h"

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
TimeAlpha::TimeAlpha(const char * name, std::string keylist ) : BCModel(name) {
	// constructor
	// define parameters here. For example:
	// AddParameter("mu",-2,1,"#mu");
	// and set priors, if using built-in priors. For example:
	// SetPriorGauss("mu",-1,0.25);
	ReadData( keylist );
}

// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
}

// ---------------------------------------------------------
// Read Data from Runs key list tier3
// ---------------------------------------------------------
int TimeAlpha::ReadData( string keylist )
{
	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");

	//
	gada::FileMap myMap;
	myMap.SetRootDir( GERDA_PHASEII_DATA );
	myMap.BuildFromListOfKeys( keylist );

	gada::DataLoader l;
	l.AddFileMap(&myMap);
	l.BuildTier3();

	TChain * chain = l.GetSharedMasterChain();
	int nentries = chain->GetEntries();

	cout << "There are " << nentries << " events in the chain!" <<endl;

	// fill the data in histograms
	int eventChannelNumber;
	unsigned long long timestamp;
	unsigned int decimalTimestamp;
	vector<int> * firedFlag = new vector<int>(fNDetectors);
	int multiplicity;
	vector<double> * energy = new vector<double>(fNDetectors);
	int isTP;
	int isVetoed;
	int isVetoedInTime;
	vector<int> * failedFlag = new vector<int>(fNDetectors);

	chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
	chain -> SetBranchAddress("timestamp",&timestamp);
	chain -> SetBranchAddress("decimalTimestamp",&decimalTimestamp);
	chain -> SetBranchAddress("firedFlag", &firedFlag);
	chain -> SetBranchAddress("multiplicity",&multiplicity);
	chain -> SetBranchAddress("rawEnergyGauss",&energy);
	chain -> SetBranchAddress("isVetoed", &isVetoed);
	chain -> SetBranchAddress("isTP",&isTP);
	chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
	chain -> SetBranchAddress("failedFlag",&failedFlag);

	chain->GetEntry( nentries - 1 );

	// Prepare histograms cut the last piece of data if it doen't fit in the 20 days scheme
	unsigned long secondIn20Days = 60*60*24*20;
	unsigned long nBins = timestamp / secondIn20Days;
	int maxDay = 20 * nBins;

	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", nBins, 0, maxDay );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", nBins, 0, (double)maxDay );

	for( int e = 0; e < nentries; e++ )
	{
		chain->GetEntry(e);

		// Apply cuts
		if( multiplicity != 1 ) continue;
		if( isVetoed ) 			continue;
		if( isVetoedInTime ) 	continue;
		if( failedFlag ) 		continue;

		if( isTP ) fHLiveTimeFraction->Fill( timestamp / secondIn20Days );
		else if( energy->at(eventChannelNumber) > 3500. && energy->at(eventChannelNumber) < 5300. ) 
			fHTimeAlpha->Fill( timestamp / secondIn20Days );
	}

	// Expected number of TP in 20 days: Test Pulser rate is 0.05Hz = 1/20s
	int TPExpected = 24 * 60 * 60 ;

	fHLiveTimeFraction->Scale( 1./(double)TPExpected );

	return 0;
}


// ---------------------------------------------------------
double TimeAlpha::LogLikelihood(const std::vector<double> & parameters) {
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	// access parameters from vector by remembering their positions, e.g.
	// double mu = parameters[0];
	// or by looking up their indicies, e.g.
	// double mu = parameters[fParameters.Index("mu")];

	// Calculate your likelihood according to your model. You may find
	// the built in functions such as BCMath::LogPoisson helpful.
	// Return the logarithm of this likelood

	return -1;
}

// ---------------------------------------------------------
// double TimeAlpha::LogAPrioriProbability(const std::vector<double> & parameters) {
// 	// This method returns the logarithm of the prior probability for the
// 	// parameters p(parameters).

// 	// You need not overload this function, if you are using built-in
// 	// priors through the function SetPriorGauss, SetPriorConstant, etc.
// }
