// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <string>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TimeAlpha.h"

#include <BAT/BCMath.h>
#include <BAT/BCH1D.h>

#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
TimeAlpha::TimeAlpha( const char * name ) : BCModel(name) {

	fDetectorType = all;
	fModel = LinAndExp;

	fNBins = 100;
	fHMaximum = 0.;
	fHMinimum = 100.;

	fHTimeAlpha = NULL;
	fHLiveTimeFraction = NULL;

	fVTimeAlpha = vector<double>(0);
	fVLiveTimeFraction = vector<double>(0);

	// constructor
	// define parameters here. For example:
	AddParameter( "constant", 0., 1000. );
	AddParameter( "amplitude", 0., 1000. );
	AddParameter( "halflife", 0., 1000. );

	// and set priors, if using built-in priors. For example:
	switch ( fModel )
	{
		case LinAndExp:
			SetPriorConstant( 0 );
			SetPriorConstant( 1 );
		//	SetPriorConstant( 2 );
			SetPriorGauss( 2, 138.3763, 40. );
			break;
	}

	cout << "Created model" << endl;
}

// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
	if(fHTimeAlpha) fHTimeAlpha->Delete();
	if(fHLiveTimeFraction) fHLiveTimeFraction->Delete();
}

// ---------------------------------------------------------
// Read Data from Runs key list tier3
// ---------------------------------------------------------
int TimeAlpha::ReadData( string keylist )
{
	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");
	string GERDA_DATA_SETS = getenv("GERDA_DATA_SETS"); GERDA_DATA_SETS += "/";
	string data_set = GERDA_DATA_SETS; data_set += keylist;

	cout << "Reading data from keylist: " << data_set << endl;
	cout << "This may take a while..." << endl;

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

	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", fNBins, fHMinimum, fHMaximum );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", fNBins, fHMinimum, fHMaximum );

	TH1D * HTimeAlpha_fine = new TH1D( "HTimeAlpha_fine", "HTimeAlpha", (int)fHMaximum, fHMinimum, fHMaximum );
	TH1D * HLiveTimeFraction_fine = new TH1D( "HLiveTimeFraction_fine", "HLiveTimeFraction", (int)fHMaximum, fHMinimum, fHMaximum );

	cout << "Loop over event chain." << endl;

	vector<int> enrBEGeChannels = {	0,1,2,3,4,5,6,7,
									11,12,13,14,15,16,17,18,
									19,20,21,22,23,24,25,26,
									30,31,32,33,34,35 };
	vector<int> enrCoaxChannels = {8,9,10,27,28,29,36};
	vector<int> natCoaxChannels = {37,38,39};

	for( int e = 0; e < nentries; e++ )
	{
		chain->GetEntry(e);

		// Apply cuts
		// if( multiplicity != 1 ) continue;
		// if( isVetoed ) 			continue;
		// if( isVetoedInTime ) 	continue;

		double time = (double)timestamp / (double)secondsInAnHour;
		time -= time0;

		if( e%10000 == 0 ) cout << "h: " << time << " d: " << time/24. << endl;


		for( int i = 0; i < fNDetectors; i++ )
		{
			if( failedFlag->at(i) ) continue;

			if( fDetectorType == enrBEGe )
			{
				if( find( enrBEGeChannels.begin(), enrBEGeChannels.end(), i) == enrBEGeChannels.end() )
					continue;
			}
			if( fDetectorType == enrCoax )
			{
				if( find( enrCoaxChannels.begin(), enrCoaxChannels.end(), i) == enrCoaxChannels.end() )
					continue;
			}
			if( fDetectorType == natCoax )
			{
				if( find( natCoaxChannels.begin(), natCoaxChannels.end(), i) == natCoaxChannels.end() )
					continue;
			}

			if( isTP )
			{
				fHLiveTimeFraction->Fill( time / 24. );
				HLiveTimeFraction_fine->Fill( time / 24. );

				break;
			}
			else if( energy->at(i) > 3500. && energy->at(i) < 5300. )
			{
				fHTimeAlpha->Fill( time / 24. );
				HTimeAlpha_fine->Fill( time / 24. );
			}
		}
	}

	// Expected number of TP in 1 day: Test Pulser rate is 0.05Hz = 1/20s
	double TPExpected = 24. * 60. * 60. / 20.;

	fHLiveTimeFraction->Scale( 1./( TPExpected * fHTimeAlpha -> GetBinWidth( 1 ) ) );
	HLiveTimeFraction_fine->Scale( 1./( TPExpected * HLiveTimeFraction_fine-> GetBinWidth( 1 ) ) );

	TFile * out = new TFile( "test.root", "RECREATE" );
	HTimeAlpha_fine -> Write();
	fHTimeAlpha -> Write();
	HLiveTimeFraction_fine -> Write();
	fHLiveTimeFraction -> Write();
	out -> Close();

	FillDataArrays();

	return 0;
}


// ---------------------------------------------------------
int TimeAlpha::FillDataArrays()
{
	cout << "Loading data in vectors for faster access." << endl;

	for( int b = 1; b <= fNBins; b++ )
	{
		fVTimeAlpha.push_back( fHTimeAlpha->GetBinContent(b) );
		fVLiveTimeFraction.push_back( fHLiveTimeFraction->GetBinContent(b) );
	}

	return 0;
}


// ---------------------------------------------------------
double TimeAlpha::LogLikelihood(const std::vector<double> & parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	// access parameters from vector by remembering their positions, e.g.
	// double mu = parameters[0];
	// or by looking up their indicies, e.g.
	// double mu = parameters[fParameters.Index("mu")];

	// Calculate your likelihood according to your model. You may find
	// the built in functions such as BCMath::LogPoisson helpful.
	// Return the logarithm of this likelood

	double constant = parameters[0];
	double amplitude = parameters[1];
	double halflife = parameters[2];

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

	double logprob = 0.;

	for( int b = 0; b < fNBins; b++)
	{
		double t1 = fHMinimum + b * BinWidth;
		double t2 = t1 + BinWidth;

		double tau = halflife / log(2.);

		double lambda = constant * BinWidth + amplitude * tau * ( exp( -t2/tau ) - exp( -t1/tau ) );
		lambda *= fVLiveTimeFraction.at(b);

	    	double data = fVTimeAlpha.at(b);

	    	double sum = data * log(lambda) - lambda - BCMath::LogFact( (int)data );

		if( lambda > 0 ) logprob += sum;
	  }

	  return logprob;
}



// ---------------------------------------------------------
void TimeAlpha::WriteOutput( string outputfilename )
{
	// Write Livetime Fraction

	// Write Data

	// Write Model
	const vector<double> BestFit = GetBestFitParametersMarginalized();

	vector<double> * BestFitParameters;
	vector<double> * BestFitParameterErrors1s;
//	vector<double> * BestFitParameterErrors2s;
//	vector<double> * BestFitParameterErrors3s;

	for( uint i = 0; i < GetNParameters(); i++ )
	{
		BCH1D * parHisto = MCMCGetH1Marginalized( i );

		double value =  parHisto->GetMode();

		BestFitParameters -> push_back( value );

		double xmin, xmax;
	    parHisto->GetSmallestInterval( xmin, xmax, 0.683 );
	    BestFitParameterErrors1s -> push_back( (xmax - xmin) / 2. );

	    cout << "Par0 " << value << "(" << xmin << "," << xmax << ")" << endl;

	    /*
	    parHisto->GetSmallestInterval( xmin, xmax, 0.954 );
	    BestFitParameterErrors2s -> push_back( (xmax - xmin) / 2. );
	    parHisto->GetSmallestInterval( xmin, xmax, 0.997 );
	    BestFitParameterErrors3s -> push_back( (xmax - xmin) / 2. );
	    */
	}

	TF1 * model = new TF1( "model", "[0] + [1]*exp( -x/[3] )", fHMinimum, fHMaximum );
	model -> SetParameters( BestFitParameters->at(0),
			BestFitParameters->at(1), BestFitParameters->at(2) );
/*	model -> SetParErrors( 0, BestFitParameterErrors1s(0) );
	model -> SetParErrors( 1, BestFitParameterErrors1s(1) );
	model -> SetParErrors( 2, BestFitParameterErrors1s(2) );
*/

	TH1D * expected = new TH1D( "histo_model", "histo_model", fNBins,  fHMinimum, fHMaximum );

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

	for( int b = 0; b < fNBins; b++)
	{
		double t1 = fHMinimum + b * BinWidth;
		double t2 = fHMinimum + (b+1) * BinWidth;

		double lambda = BestFitParameters->at(0) * ( t2 - t1 )
				- BestFitParameters->at(1) / BestFitParameters->at(2)
				* ( exp( -t2/BestFitParameters->at(2) ) - exp( -t1/BestFitParameters->at(2) ) );
		lambda *= fVLiveTimeFraction.at(b);

		expected -> SetBinContent( b, lambda );
	}

	TFile * out = new TFile( outputfilename.c_str(), "RECREATE" );
	model -> Write();
	fHTimeAlpha -> Write();
	fHLiveTimeFraction -> Write();
	expected -> Write();
	out -> Close();

	return;
}


// ---------------------------------------------------------
void TimeAlpha::WriteDistributions()
{	
	TFile * test = new TFile( "posteriors.root", "RECREATE" );

	for( uint i = 0; i < GetNParameters(); i++ )
	{
		BCH1D * parHisto = MCMCGetH1Marginalized( i );

		TH1D * histo = parHisto -> GetHistogram();

		histo -> Write();
	}
	
	test -> Close();

	return;
}


// ---------------------------------------------------------
// double TimeAlpha::LogAPrioriProbability(const std::vector<double> & parameters) {
// 	// This method returns the logarithm of the prior probability for the
// 	// parameters p(parameters).

// 	// You need not overload this function, if you are using built-in
// 	// priors through the function SetPriorGauss, SetPriorConstant, etc.
// }
