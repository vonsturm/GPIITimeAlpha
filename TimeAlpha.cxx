// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <string>
#include <climits>

#include "TimeAlpha.h"

#include <BAT/BCMath.h>

#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"

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
	string GERDA_DATA_SETS = getenv("GERDA_DATA_SETS"); GERDA_DATA_SETS += "/";
	string data_set = GERDA_DATA_SETS; data_set += keylist;

	//
	gada::FileMap myMap;
	myMap.SetRootDir( GERDA_PHASEII_DATA );
	myMap.BuildFromListOfKeys( data_set );

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

	// Prepare histograms cut the last piece of data if it doen't fit in the 20 days scheme
	unsigned long long secondsInAnHour= 60 * 60;

	chain->GetEntry(0);

	double time0 = (double)timestamp / (double)secondsInAnHour; // time in hours

	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", 10, 0, 200. );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", 10, 0, (double)200. );
	fHRealDecay = new TH1D( "HRealDecay", "HRealDecay", 10, 0, (double)200. );

	for( int e = 0; e < nentries; e++ )
	{
		chain->GetEntry(e);

		// Apply cuts
		// if( multiplicity != 1 ) continue;
		if( isVetoed ) 			continue;
		if( isVetoedInTime ) 	continue;

		double time = (double)timestamp / (double)secondsInAnHour;
		time -= time0;

		if( e%10000 == 0 ) cout << "h: " << time << " d: " << time/24. << endl;


		for( int i = 0; i < fNDetectors; i++ )
		{
			if( failedFlag->at(i) ) continue;

			if( isTP ) fHLiveTimeFraction->Fill( time / 24. );
			else if( energy->at(i) > 3500. && energy->at(i) < 5300. )
				fHTimeAlpha->Fill( time / 24. );

			break;
		}
	}

	// Expected number of TP in 20 days: Test Pulser rate is 0.05Hz = 1/20s
	int TPExpected = 24 * 60 * 60 ;

	fHLiveTimeFraction->Scale( 1./(double)TPExpected );

	fHRealDecay->Add( fHTimeAlpha );
	fHRealDecay->Divide( fHLiveTimeFraction );

	TFile * rootfile = new TFile( "test.root", "RECREATE" );
	fHTimeAlpha->Write();
	fHLiveTimeFraction->Write();
	fHRealDecay->Write();
	rootfile->Close();

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
