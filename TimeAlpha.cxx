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
#include "TCanvas.h"

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

#include "GETRunConfiguration.hh"
#include "GERunConfigurationManager.hh"

#include "ProgressBar.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
TimeAlpha::TimeAlpha( const char * name ) : BCModel(name) {

	fDetectorType = kUNKNOWN; // all detectors
	fDetectorEnriched = true;
//	fModel = OnlyLin;
//	fModel = OnlyExp;
	fModel = LinAndExp;

	fNBins = 100;
	fHMaximum = 0.;
	fHMinimum = 100.;

	fHTimeAlpha = NULL;
	fHLiveTimeFraction = NULL;

	fVTimeAlpha = vector<double>(0);
	fVLiveTimeFraction = vector<double>(0);

	DefineParametersAndPriors();

	cout << "Created model" << endl;
}

int TimeAlpha::DefineParametersAndPriors()
{
	// constructor
	// define parameters here. For example:
	AddParameter( "constant", 0., 15. );
	AddParameter( "amplitude", 0., 25. );
	AddParameter( "halflife", 137., 140. );

	// and set priors, if using built-in priors. For example:
	switch ( fModel )
	{
		case LinAndExp:
			fModelName = "linexp";
			SetPriorConstant( 0 );
			SetPriorConstant( 1 );
//			SetPriorConstant( 2 );
			SetPriorGauss( 2, 138.4, 0.2 );
//			SetPriorDelta( 2, 138.3763 );
			break;

		case OnlyExp:
			fModelName = "exp";
			SetPriorDelta( 0, 0. );
			SetPriorConstant( 1 );
			SetPriorConstant( 2 );
//			SetPriorGauss( 2, 138.4, 0.2 );
			break;

		case OnlyLin:
			fModelName = "lin";
			SetPriorConstant( 0 );
			SetPriorDelta( 1, 0. );
			SetPriorDelta( 2, 0. );
			break;
		default:
			cout << "Model not set!" << endl;
			return -1;
	}

	return 0;
}


// ---------------------------------------------------------
TimeAlpha::~TimeAlpha() {
	// destructor
	if(fHTimeAlpha) fHTimeAlpha->Delete();
	if(fHLiveTimeFraction) fHLiveTimeFraction->Delete();
}

// ---------------------------------------------------------
int TimeAlpha::InitializeHistograms()
{
	ResetHistograms();

	fHTimeAlpha = new TH1D( "HTimeAlpha", "HTimeAlpha", fNBins, fHMinimum, fHMaximum );
	fHLiveTimeFraction = new TH1D( "HLiveTimeFraction", "HLiveTimeFraction", fNBins, fHMinimum, fHMaximum );

	fHTimeAlpha_fine = new TH1D( "HTimeAlpha_fine", "HTimeAlpha_fine", (int)fHMaximum, fHMinimum, fHMaximum );
	fHLiveTimeFraction_fine = new TH1D( "HLiveTimeFraction_fine", "HLiveTimeFraction_fine", (int)fHMaximum, fHMinimum, fHMaximum );

	if( fHTimeAlpha && fHLiveTimeFraction && fHTimeAlpha_fine && fHLiveTimeFraction_fine )
		return 0;
	else
		return -1;
}

// ---------------------------------------------------------
// Read Data from Runs key list tier3
// ---------------------------------------------------------
int TimeAlpha::ReadDataPhaseII( string keylist )
{
	string MU_CAL = getenv("MU_CAL");

	// check if $MU_CAL is set
	if( MU_CAL == "" )
	{
		cout << "Environment variable MU_CAL not set. Have it point to the directory where runconfiguration.db can be found.\n"
			 << "On LNGS that would be e.g. /nfs/gerda5/gerda-data/blind/active/meta/config/_aux/geruncfg " << endl;
		return -1;
	}

	cout << "MU_CAL = " << MU_CAL << endl;

	// initialize run configuration manager
	GERunConfigurationManager * RunConfManager = new GERunConfigurationManager();
	GETRunConfiguration * RunConf = 0;
	RunConfManager -> AllowRunConfigurationSwitch(true);
	RunConfManager -> SetVerbosity(1);

	// tell the data map where to find everything
	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");
	GERDA_PHASEII_DATA += "/gen";
	string GERDA_DATA_SETS = getenv("GERDA_DATA_SETS"); GERDA_DATA_SETS += "/";
	string data_set = GERDA_DATA_SETS; data_set += keylist;

	cout << "Reading from data dir: " << GERDA_PHASEII_DATA << endl;
	cout << "\t data set dir: " << GERDA_DATA_SETS << endl;
	cout << "\t from keylist: " << keylist << endl;
	cout << "This may take a while..." << endl;

	// setting FileMap and DataLoader
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
	//unsigned int decimalTimestamp;
	int multiplicity;
	int isTP;
	int isVetoedInTime;
	vector<int> * firedFlag = new vector<int>(fNDetectors);
	vector<double> * energy = new vector<double>(fNDetectors);
	vector<double> * risetime = new vector<double>(fNDetectors);
	vector<int> * failedFlag = new vector<int>(fNDetectors);
	vector<int> * failedFlag_isPhysical = new vector<int>(fNDetectors);
	vector<int> * failedFlag_isSaturated = new vector<int>(fNDetectors);

	chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
	chain -> SetBranchAddress("timestamp",&timestamp);
	//chain -> SetBranchAddress("decimalTimestamp",&decimalTimestamp);
	chain -> SetBranchAddress("multiplicity",&multiplicity);
	chain -> SetBranchAddress("isTP",&isTP);
	chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
	chain -> SetBranchAddress("firedFlag", &firedFlag);
	chain -> SetBranchAddress("rawEnergyGauss",&energy);
	chain -> SetBranchAddress("risetime",&risetime);
	chain -> SetBranchAddress("failedFlag_isPhysical",&failedFlag_isPhysical);
	chain -> SetBranchAddress("failedFlag_isSaturated",&failedFlag_isSaturated);

	chain -> GetEntry( 0 ); // Get timestamp of first event in chain
	unsigned long long time0 = timestamp; // time in hours

	chain -> GetEntry( nentries ); // Get timestamp of last event in chain
	unsigned long long timeN = timestamp; // time in hours

	// calculate fit interval in days
	int timeInDays = (int)( ( timeN - time0 ) / ( 60. * 60. * 24. ) );

	SetNBinsHistograms( timeInDays/fBinning + 1, 0., timeInDays - timeInDays%fBinning + fBinning );
	InitializeHistograms();

	cout << "Fitting " << timeInDays << "days of data." << endl;

	cout << "Starting loop over event chain..." << endl;

	// ProgressBar
	ProgressBar bar( nentries, '#', false );

	for( int e = 0; e < nentries; e++ )
	{
		bar.Update();

		chain->GetEntry(e);

		unsigned long long time = timestamp - time0;

		if( isTP )
		{
			fHLiveTimeFraction -> Fill( time );
			fHLiveTimeFraction_fine -> Fill( time );
			continue;
		}

		// Apply cuts
		if( isVetoedInTime ) continue;
		if( multiplicity > 1 ) continue;

		for( int d = 0; d < fNDetectors; d++ )
		{
			if( !firedFlag -> at( d ) ) continue;
			if( failedFlag -> at( d ) ) continue;

			// use kUNKNOWN to select all detectors
			RunConf = RunConfManager -> GetRunConfiguration( timestamp );
			DetectorType_t dType = RunConf -> GetDetector( d ) -> GetDetectorType(); // kIsBEGe=1, kIsCOAX=2, kUNKNOWN=3
			Bool_t dIsEnriched = RunConf -> GetDetector( d ) -> IsEnriched();
			string dDetName = RunConf -> GetDetector( d ) -> GetDetName();

			if( !RunConf -> IsOn( d ) ) continue;
			if( fSingleDetectorFit ){ if( dDetName != fDataSet ) continue; }
			else if( fDetectorType != kUNKNOWN && fDetectorEnriched != kUNKNOWN )
			{
				if( dType != fDetectorType ) continue;
				if( dIsEnriched != fDetectorEnriched ) continue;
			}

			if( multiplicity == 1 && !failedFlag_isPhysical -> at( d ) && energy -> at( d ) > 3500.
				|| !failedFlag_isSaturated -> at( d ) )
			{
				fHTimeAlpha -> Fill( time );
				fHTimeAlpha_fine -> Fill( time );

				if(  failedFlag_isSaturated -> at( d ) )
					cout << "Saturated Event: " << e << "(" << d << ")" << endl;
			}
		}
	}

	cout << "...DONE" << endl;

	// Test Pulser rate is 0.05Hz = 1/20s
	double TPFrequency = 1. / 20.;

	// Scale live time fraction with TP rate
	fHLiveTimeFraction -> Scale( 1./( TPFrequency * fHLiveTimeFraction -> GetBinWidth( 1 ) ) );
	fHLiveTimeFraction_fine -> Scale( 1./( TPFrequency * fHLiveTimeFraction_fine -> GetBinWidth( 1 ) ) );

	TFile * out = new TFile( "./out/TimeAlpha_Data.root", "RECREATE" );
	fHTimeAlpha -> Write();
	fHTimeAlpha_fine -> Write();
	fHLiveTimeFraction -> Write();
	fHLiveTimeFraction_fine -> Write();
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
		if( fVLiveTimeFraction.at(b) <= fLTLimit ) continue;

		double t1 = fHMinimum + b * BinWidth;
		double t2 = t1 + BinWidth;

		double tau = halflife / log(2.);

		double lambda = constant * BinWidth - amplitude * tau * ( exp( -t2/tau ) - exp( -t1/tau ) );
		lambda *= fVLiveTimeFraction.at(b);

		if( lambda <= 0. ) continue;

	    double data = fVTimeAlpha.at(b);

	    double sum = data * log(lambda) - lambda - BCMath::LogFact( (int)data );

		logprob += sum;
	  }

	  return logprob;
}


// ---------------------------------------------------------
double TimeAlpha::EstimatePValue()
{
	//Allen's routine for the evaluation of p-value
  	//This is derived from PRD 83 (2011) 012004, appendix
  	//taken from Luciano

	double logp0 = LogLikelihood( GetBestFitParameters() ); //likelihood at the mode

  	/*
    	Now we initialize the Markov Chain by setting the number of entries in a bin to
    	the integer
    	part of the mean in that bin.  This will give the maximum possible likelihood.
  	*/

	double sumlog = 0;

	/* mean is the array where you have the expected mean in each bin (calculated
     	from the parameter values at the mode. Nom is the nominal value of entries in
     	each bin (integer part)
  	*/

  	vector<double> mean( fNBins, 0 );
	vector<int> nom( fNBins, 0 );

  	vector<double> parameters = GetBestFitParameters();
	double constant = parameters[0];
	double amplitude = parameters[1];
	double halflife = parameters[2];

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

  	for( int ibin = 0; ibin < fNBins; ibin++ )
	{
		//Require a live time fraction of at least 10% otherwise skip data point
		if( fVLiveTimeFraction.at(ibin) <= fLTLimit ) continue;

                double t1 = fHMinimum + ibin * BinWidth;
                double t2 = t1 + BinWidth;

                double tau = halflife / log(2.);

                double lambda = constant * BinWidth - amplitude * tau * ( exp( -t2/tau ) - exp( -t1/tau ) );
                lambda *= fVLiveTimeFraction.at(ibin);

		mean.at( ibin ) = std::max( lambda, 1e-8 );
      		nom.at( ibin ) = int( mean.at( ibin ) );
      		sumlog += BCMath::LogPoisson( nom[ ibin ], mean[ ibin ] );
  	}

  	cout << "Logprob for best: " << sumlog << endl;

  	/*
  	Now we run the Markov chain to generate new configurations.  One iteration means
  	a loop over all the bins, with an attempt to vary each bin up or down by one unit.  We
  	accept/reject at each step  and compare the data logprob to the simulated at the end of each iteration.
  	*/

  	const int nloops = 100000;
  	int Pgood = 0;

  	for( int iloop = 0; iloop < nloops; iloop++ )
  	{
      		for( int ibin = 0; ibin < fNBins; ibin++)
      		{
			if( fVLiveTimeFraction.at(ibin) <= fLTLimit ) continue;

    	  		if ( rand() > RAND_MAX/2 ) // Try to increase the bin content by 1
    	  		{
    		  		double r = mean[ ibin ]/( nom[ ibin ] + 1 );
    		  		double rtest = double( rand() )/RAND_MAX;
    		  		if( rtest < r ) //Accept
    		  		{
    			  		nom[ ibin ] = nom[ ibin ] + 1;
    			  		sumlog += log(r);
    		  		}
    	  		}
    	  		else // Try to decrease the bin content by 1
    	  		{
    		  		double r = nom[ ibin ]/mean[ ibin ];
    		  		double rtest = double( rand() )/RAND_MAX;
    		  		if ( rtest < r ) //Accept
    		  		{
    			  		nom[ ibin ] = nom[ ibin ] - 1;
    			  		sumlog += log(r);
    		  		}
    	  		}
      		}

		if ( sumlog < logp0 ) Pgood++;
  	}

  	double pvalue = double(Pgood)/double(nloops);

  	cout << "p-value is " << pvalue << endl;

  	return pvalue;
}



// ---------------------------------------------------------
void TimeAlpha::WriteOutput( string outputfilename )
{
	// Write Livetime Fraction

	// Write Data

	// Write Model
	const vector<double> BestFit = GetBestFitParametersMarginalized();

	vector<double> * BestFitParameters = new vector<double>;
	vector<double> * BestFitParameterErrors1s = new vector<double>;
//	vector<double> * BestFitParameterErrors2s = new vector<double>;
//	vector<double> * BestFitParameterErrors3s = new vector<double>;

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

	TF1 * modelfunction = new TF1( "modelfunction", "[0]+[1]*exp(-x*log(2.)/[2])", fHMinimum, fHMaximum );
	modelfunction -> SetParameters( BestFitParameters->at(0), BestFitParameters->at(1), BestFitParameters->at(2) );
	modelfunction -> SetParError( 0, BestFitParameterErrors1s->at(0) );
	modelfunction -> SetParError( 1, BestFitParameterErrors1s->at(1) );
	modelfunction -> SetParError( 2, BestFitParameterErrors1s->at(2) );

	cout << fHMinimum << " " << fHMaximum << endl;
	cout << BestFitParameters->at(0) << " " <<	BestFitParameters->at(1) << " "<<	BestFitParameters->at(2) << " " << endl;
	cout << BestFitParameterErrors1s->at(0) << " " <<	BestFitParameterErrors1s->at(1) << " "<<	BestFitParameterErrors1s->at(2) << " " << endl;

	TCanvas * c = new TCanvas("c1");
	modelfunction -> Draw();

	TH1D * expected = new TH1D( "histo_model", "histo_model", fNBins,  fHMinimum, fHMaximum );

	double BinWidth = (fHMaximum - fHMinimum) / fNBins;

	for( int b = 0; b < fNBins; b++)
	{
		double t1 = fHMinimum + b * BinWidth;
		double t2 = t1 + BinWidth;

		double tau = BestFitParameters->at(2) / log(2.);

		double lambda = BestFitParameters->at(0) * BinWidth
				- BestFitParameters->at(1) * tau
				* ( exp( -t2/tau ) - exp( -t1/tau ) );
		lambda *= fVLiveTimeFraction.at(b);

		expected -> SetBinContent( b+1, lambda );
	}

	TFile * out = new TFile( outputfilename.c_str(), "RECREATE" );
	c -> Write();
	modelfunction -> Write();
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
